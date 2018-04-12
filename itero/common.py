
import os
import sys
import shutil
import tarfile
import numpy

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from itero import samtools
from itero import bedtools
from itero import spades

import pdb


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
        except AssertionError:
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
            assembly_fasta_fname = os.path.join(sample_dir_iter, "loci","{}.fasta".format(locus))
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
    spades_assembly_dir = spades.spades_paired_end_assembly(iteration, sample, sample_dir_iter_locus, fastqs, locus, clean)
    spades_assembly_fasta = os.path.join(spades_assembly_dir, "contigs.fasta")
    # if the assembly exists, copy it
    if os.path.isfile(spades_assembly_fasta):
        shutil.copyfile(spades_assembly_fasta, os.path.join(sample_dir_iter, "loci","{}.fasta".format(locus)))
    if clean:
        shutil.rmtree(sample_dir_iter_locus)
    sys.stdout.write('.')
    sys.stdout.flush()
    return 0


# this works ok, except we need to only be zipping the individual locus files.
def zip_assembly_dir(log, sample_dir_iter, clean, prev_iter):
    prev_sample_locus_iter = os.path.join(get_previous_sample_dir_iter(log, sample_dir_iter, prev_iter), "loci")
    if not clean:
        log.info("Zipping the locus directory for iter-{}. May be slow, particularly on HPC.".format(prev_iter))
        #pdb.set_trace()
        output_tarfile = "{}.tar.gz".format(prev_sample_locus_iter)
        with tarfile.open(output_tarfile, "w:gz") as tar:
            tar.add(prev_sample_locus_iter, arcname=os.path.basename(prev_sample_locus_iter))
        # remove unzipped directory
        shutil.rmtree(prev_sample_locus_iter)
    else:
        # NUKE IT
        # remove unzipped directory
        shutil.rmtree(prev_sample_locus_iter)