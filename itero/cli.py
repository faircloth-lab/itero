#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2017 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 17 July 2017 13:35 CDT (-0500)
"""


import os
import sys
import time
import glob
import shutil
import tarfile
import argparse
import subprocess
import ConfigParser
import multiprocessing

import numpy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mpi4py import MPI

from itero import bwa
from itero import samtools
from itero import bedtools
from itero import spades

from itero.helpers import FullPaths, CreateDir, is_dir, is_file
from itero.raw_reads import get_input_files
from itero.log import setup_logging

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Iteratively assemble loci from raw reads and a seed file"""
    )
    parser.add_argument(
        "--config",
        type=is_file,
        action=FullPaths,
        default=None,
        help="""A configuration file containing reads to assemble"""
    )
    parser.add_argument(
        "--subfolder",
        type=str,
        default='',
        help="""A subdirectory, below the level of the group, containing the reads"""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=CreateDir,
        help="""The directory in which to store the output"""
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=5,
        help="""The number of iterations to run for each locus"""
    )
    parser.add_argument(
        "--local-cores",
        type=int,
        default=1,
        help="""The number of cores to use on the main node"""
    )
    parser.add_argument(
        "--use-mpi",
        action="store_true",
        default=False,
        help="""Perform assemblies using MPI.  Number of cores passed with mpirun.""",
    )
    parser.add_argument(
        "--mpi-cores",
        type=int,
        default=None,
        help="""The number of cores to use for MPI"""
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        default=False,
        help="""Cleanup all intermediate files""",
    )
    parser.add_argument(
        "--only-single-locus",
        action="store_true",
        default=False,
        help="""Assemble only to a single contig""",
    )
    parser.add_argument(
        "--allow-multiple-contigs",
        action="store_true",
        default=False,
        help="""Allow assembly stages to produce multiple contigs""",
    )
    parser.add_argument(
        "--do-not-zip",
        action="store_true",
        default=False,
        help="""Do not zip the iteration files, which is the default behavior.""",
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="""The logging level to use."""
    )
    parser.add_argument(
        "--log-path",
        action=FullPaths,
        type=is_dir,
        default=None,
        help="""The path to a directory to hold logs."""
    )
    p = parser.parse_args()
    if p.use_mpi and p.mpi_cores is None:
        parser.error('You specified --use-mpi which means you must also specify --mpi-cores <int>')
    else:
        return p


def get_input_data(log, args, conf):
    individuals = conf.items('individuals')
    updated_individuals = []
    for sample in individuals:
        try:
            # deal with relative paths in config
            if sample[1].startswith(".."):
                pth = os.path.join(os.path.dirname(args.config), sample[1])
            else:
                pth = sample[1]
            assert os.path.isdir(pth)
            updated_individuals.append((sample[0], pth))
        except:
            raise IOError("{} is not a directory".format(sample[1]))
    #return reference, individuals
    return updated_individuals


def get_seed_names(seeds):
    with open(seeds, "ru") as infile:
        return [i.lstrip(">").rstrip() for i in infile if i.startswith(">")]



def get_fasta(log, sample, sample_dir_iter, locus_names, multiple_hits=False, iteration=0):
    assemblies = []
    assemblies_stats = []
    all_fasta_out_fname = os.path.join(sample_dir_iter, 'iter-{}.all-fasta.fasta'.format(iteration))
    all_fasta_stats_fname = os.path.join(sample_dir_iter, 'iter-{}.all-fasta.stats.csv'.format(iteration))
    print("")
    for locus in locus_names:
        try:
            assembly_fasta_fname = os.path.join(sample_dir_iter, "loci", locus, "{}-assembly".format(locus), "contigs.fasta")
            sequence = list(SeqIO.parse(assembly_fasta_fname, 'fasta'))
            if len(sequence) > 1 and multiple_hits is True:
                # keep only contigs > 100 bp
                for v, seq in enumerate(sequence):
                    #pdb.set_trace()
                    if len(seq) >= 100:
                        if v == 0:
                            new_seq = SeqRecord(seq.seq)
                            new_seq.id = seq.id.replace("NODE", locus.split("_")[0])
                            new_seq.description = ""
                            new_seq.name = ""
                        else:
                            #pdb.set_trace()
                            new_seq.seq += Seq("{}".format(200*'N')) + seq.seq
                    else:
                        pass
                log.warn("Locus {} has multiple hits (allowed during initial rounds).  Padded both with Ns and put back into seeds.".format(new_seq.id.split("_")[0]))
                assemblies.append(new_seq)
                assemblies_stats.append(len(new_seq.seq.strip("N")))
            elif (len(sequence) == 1 or multiple_hits is False) and len(sequence[0]) >= 100:
                seq = sequence[0]
                seq.id = seq.id.replace("NODE", locus.split("_")[0])
                seq.description = ""
                seq.name = ""
                assemblies.append(seq)
                assemblies_stats.append(len(seq))
            else:
                log.warn("Dropped locus {} for having multiple contigs or being short (<100 bp)".format(locus))
        except IOError:
            log.warn("Dropped locus {} for having no assembled contigs (or coverage < 5)".format(locus))
    if len(assemblies) == 0:
        log.critical("Zero valid contigs were assembled.  Quitting.")
        sys.exit()
    with open(all_fasta_out_fname, 'w') as outfile:
        for assembly in assemblies:
            try:
                outfile.write(assembly.format('fasta'))
            except TypeError:
                log.error("Dropped sequence {} because of unspecific sequence error".format(assembly.id.split("_")[0]))
    log.info("{} sequences. Mean sequence length {}, min {}, max {}".format(
        len(assemblies_stats),
        numpy.mean(assemblies_stats),
        numpy.min(assemblies_stats),
        numpy.max(assemblies_stats)
    ))
    numpy.savetxt(all_fasta_stats_fname, assemblies_stats, delimiter=",", fmt='%s', header='iter-{}'.format(iteration))
    return all_fasta_out_fname


def get_deltas(log, sample, sample_dir_iter, iterations, iteration=0):
    # current round of assembly
    current_fasta_stats_fname = os.path.join(sample_dir_iter, 'iter-{}.all-fasta.stats.csv'.format(iteration))
    current = numpy.genfromtxt(current_fasta_stats_fname, delimiter=",", skip_header=1)
    basename, directory = os.path.split(sample_dir_iter)
    prev_iter = get_previous_iter(log, sample_dir_iter, iterations, iteration)
    previous_sample_dir_iter = os.path.join(basename, "iter-{}".format(prev_iter))
    previous_fasta_stats_fname = os.path.join(previous_sample_dir_iter, 'iter-{}.all-fasta.stats.csv'.format(prev_iter))
    previous = numpy.genfromtxt(previous_fasta_stats_fname, delimiter=",", skip_header=1)
    difference = ((numpy.mean(current) - numpy.mean(previous)) / numpy.mean(current)) * 100
    log.info("Mean assembly length improved by {0:.0f}%".format(difference))
    return difference


def get_previous_iter(log, sample_dir_iter, iterations, iteration):
    basename, directory = os.path.split(sample_dir_iter)
    if iteration == 'final':
        # previous round of assembly
        prev_iter = iterations[-2]
    else:
        # previous round of assembly
        prev_iter = int(directory.split('-')[1]) - 1
    return prev_iter


def get_previous_sample_dir_iter(log, sample_dir_iter, prev_iter):
    basename, directory = os.path.split(sample_dir_iter)
    return os.path.join(basename, "iter-{}".format(prev_iter))


def initial_assembly(work):
    iteration, sample, sample_dir_iter, sorted_reduced_bam, locus, clean, only_single_locus = work
    sample_dir_iter_locus = os.path.join(sample_dir_iter, "loci", locus)
    os.makedirs(sample_dir_iter_locus)
    bam_paired, bam_singleton = samtools.samtools_split_bam(sample, sample_dir_iter_locus, sorted_reduced_bam, locus, clean, only_single_locus)
    fastqs = bedtools.bedtools_to_fastq(sample, sample_dir_iter_locus, bam_paired, bam_singleton, locus, clean)
    spades.spades_paired_end_assembly(iteration, sample, sample_dir_iter_locus, fastqs, locus, clean)
    sys.stdout.write('.')
    sys.stdout.flush()


# this works ok, except we need to only be zipping the individual locus files.
def zip_assembly_dir(log, sample_dir_iter, clean, prev_iter):
    log.info("Zipping the locus directory for iter-{}".format(prev_iter))
    if not clean:
        log.warn("You are not using --clean.  Zipping may be slow.")
    prev_sample_locus_iter = os.path.join(get_previous_sample_dir_iter(log, sample_dir_iter, prev_iter), "loci")
    #pdb.set_trace()
    output_tarfile = "{}.tar.gz".format(prev_sample_locus_iter)
    with tarfile.open(output_tarfile, "w:gz") as tar:
        tar.add(prev_sample_locus_iter, arcname=os.path.basename(prev_sample_locus_iter))
    # remove unzipped directory
    shutil.rmtree(prev_sample_locus_iter)

def split(container, count):
    """
    Simple function splitting a container into equal length chunks.
    Order is not preserved but this is potentially an advantage depending on
    the use case.
    """
    return [container[_i::count] for _i in range(count)]


def main():
    start_time = time.time()
    # get args and options
    args = get_args()
    # setup logging
    log, my_name = setup_logging(args)
    # get seeds from config file
    conf = ConfigParser.ConfigParser(allow_no_value=True)
    conf.optionxform = str
    conf.read(args.config)
    # get the seed file info
    seeds = conf.items("reference")[0][0]
    # deal with relative paths in config
    if seeds.startswith(".."):
        seeds = os.path.join(os.path.dirname(args.config), seeds)
    # get name of all loci in seeds file - only need to do this once
    seed_names = get_seed_names(seeds)
    # get the input data
    log.info("Getting input filenames and creating output directories")
    individuals = get_input_data(log, args, conf)
    for individual in individuals:
        sample, dir = individual
        # pretty print taxon status
        text = " Processing {} ".format(sample)
        log.info(text.center(65, "-"))
        # make a directory for sample-specific assemblies
        sample_dir = os.path.join(args.output, sample)
        os.makedirs(sample_dir)
        # determine how many files we're dealing with
        fastq = get_input_files(dir, args.subfolder, log)
        #pdb.set_trace()
        iterations = list(xrange(args.iterations)) + ['final']
        next_to_last_iter = iterations[-2]
        for iteration in iterations:
            text = " Iteration {} ".format(iteration)
            log.info(text.center(45, "-"))
            # One the last few iterations, set some things up differently to deal w/ dupe contigs.
            ## First, we'll allow multiple contigs during all but the last few rounds of contig assembly.
            ## This is because we could be assembling different parts of a locus that simply have not
            ## merged in the middle yet (but will).  We'll turn option to remove multiple contigs
            ## back on for last three rounds
            if iteration in iterations[-3:]:
                if args.allow_multiple_contigs is True:
                    allow_multiple_contigs = True
                else:
                    allow_multiple_contigs = False
            else:
                allow_multiple_contigs = True
            sample_dir_iter = os.path.join(sample_dir, "iter-{}".format(iteration))
            os.makedirs(sample_dir_iter)
            # change to sample_dir_iter
            os.chdir(sample_dir_iter)
            # copy seeds file
            if iteration == 0 and os.path.dirname(seeds) != os.getcwd():
                shutil.copy(seeds, os.getcwd())
                seeds = os.path.join(os.getcwd(), os.path.basename(seeds))
            elif iteration >= 1:
                shutil.copy(new_seeds, os.getcwd())
                seeds = os.path.join(os.getcwd(), os.path.basename(new_seeds))
            # if we are finished with it, zip the previous iteration
            if not args.do_not_zip and iteration >= 1:
                # after assembling all loci, zip the iter-#/loci directory; this will be slow if --clean is not turned on.
                prev_iter = get_previous_iter(log, sample_dir_iter, iterations, iteration)
                zipped = zip_assembly_dir(log, sample_dir_iter, args.clean, prev_iter)
            #index the seed file
            bwa.bwa_index_seeds(seeds, log)
            # map initial reads to seeds
            bam = bwa.bwa_mem_pe_align(log, sample, sample_dir_iter, seeds, args.local_cores, fastq.r1, fastq.r2, iteration)
            # reduce bam to mapping reads
            reduced_bam = samtools.samtools_reduce(log, sample, sample_dir_iter, bam, iteration=iteration)
            # remove the un-reduced BAM
            os.remove(bam)
            # sort and index bam
            sorted_reduced_bam = samtools.samtools_sort(log, sample, sample_dir_iter, reduced_bam, iteration=iteration)
            samtools.samtools_index(log, sample, sample_dir_iter, sorted_reduced_bam, iteration=iteration)
            # remove the un-sorted BAM
            os.remove(reduced_bam)
            # if we are not on our last iteration, assembly as usual
            if iteration is not 'final':
                if args.only_single_locus:
                    locus_names = ['locus-1']
                else:
                    # get list of loci in sorted bam
                    locus_names = samtools.samtools_get_locus_names_from_bam(log, sorted_reduced_bam, iteration)
                log.info("Splitting BAM and assembling")
                work = [(iteration, sample, sample_dir_iter, sorted_reduced_bam, locus_name, args.clean, args.only_single_locus) for locus_name in locus_names]
                if args.use_mpi:
                    COMM = MPI.COMM_WORLD
                    if COMM.rank == 0:
                        jobs = split(work, COMM.size)
                    else:
                        jobs = None
                    # Scatter jobs across cores.
                    jobs = COMM.scatter(jobs, root=0)
                    for job in jobs:
                        initial_assembly(job)
                else:
                    if not args.only_single_locus and args.local_cores > 1:
                        assert args.local_cores <= multiprocessing.cpu_count(), "You've specified more cores than you have"
                        pool = multiprocessing.Pool(args.local_cores)
                        pool.map(initial_assembly, work)
                    elif args.only_single_locus:
                        map(initial_assembly, work)
                    else:
                        map(initial_assembly, work)
                # after assembling all loci, get them into a single file
                new_seeds = get_fasta(log, sample, sample_dir_iter, locus_names, allow_multiple_contigs, iteration=iteration)
                # after assembling all loci, report on deltas of the assembly length
                if iteration is not 0:
                    assembly_delta = get_deltas(log, sample, sample_dir_iter, iterations, iteration=iteration)
                #
                #if iteration is 'final':
                #    prev_iter = get_previous_iter(log, sample_dir_iter, iterations, iteration)
                #    # after assembling all loci, zip the iter-#/loci directory; this will be slow if --clean is not turned on.
                #    zipped = zip_assembly_dir(log, sample_dir_iter, args.clean, prev_iter)
            elif iteration is 'final':
                log.info("Final assemblies and a BAM file with alignments to those assemblies are in {}/{}".format(individual, iteration))

    end_time = time.time()
    time_delta_sec = round(end_time - start_time, 1)
    time_delta_min = round(time_delta_sec / 60.0, 1)
    text = " Completed in {} minutes ({} seconds) ".format(time_delta_min, time_delta_sec)
    log.info(text.center(65, "="))



if __name__ == '__main__':
    main()
