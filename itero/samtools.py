#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2018 Brant Faircloth || http://faircloth-lab.org/

All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 14 April 2018 16:13 CDT (-0500)
"""

import os
import time
import subprocess

from itero.pth import get_user_path

import pdb


def samtools_version():
    cmd = [get_user_path("executables", "samtools"), '--version']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = proc.communicate()
    return stdout.split("\n")[0].split(' ')[1]


def samtools_index(log, sample, sample_dir, bam, iteration=0):
    log.info("Indexing BAM for {}".format(sample))
    cmd = [
        get_user_path("executables", "samtools"),
        "index",
        bam
    ]
    samtools_out_fname = os.path.join(sample_dir, 'iter-{}.samtools-idx.log'.format(sample))
    with open(samtools_out_fname, 'w') as samtools_out:
        proc = subprocess.Popen(cmd, stdout=samtools_out, stderr=subprocess.STDOUT)
        proc.communicate()


def samtools_reduce(log, sample, sample_dir, bam, iteration=0):
    #pdb.set_trace()
    log.info("Reducing BAM for {}, iteration {}".format(sample, iteration))
    bam_out_fname = os.path.join(sample_dir, 'iter-{}.reduce.bam'.format(iteration))
    cmd = [
        get_user_path("executables", "samtools"),
        "view",
        "-F",
        "4",
        "-bq",
        "1",
        bam,
        "-o",
        bam_out_fname
    ]
    samtools_out_fname = os.path.join(sample_dir, 'iter-{}.reduce.log'.format(iteration))
    with open(samtools_out_fname, 'w') as samtools_out:
        proc = subprocess.Popen(cmd, stdout=samtools_out, stderr=subprocess.STDOUT)
        proc.communicate()
    return bam_out_fname


def samtools_sort(log, sample, sample_dir, bam, iteration=0):
    #pdb.set_trace()
    bam_out_fname = os.path.join(sample_dir, 'iter-{}.reduce.sorted.bam'.format(iteration))
    cmd1 = [
        get_user_path("executables", "samtools"),
        "sort",
        bam,
        "-o",
        bam_out_fname
    ]
    samtools_out_fname = os.path.join(sample_dir, 'iter-{}.sort.log'.format(iteration))
    with open(samtools_out_fname, 'w') as samtools_out:
        proc = subprocess.Popen(cmd1, stdout=samtools_out, stderr=subprocess.STDOUT)
        proc.communicate()
    return bam_out_fname


def get_bam_header(log, bam, iteration):
    cmd1 = [
    get_user_path("executables", "samtools"),
    "view",
    "-H",
    bam
    ]
    proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    stdout = proc1.communicate()
    log.info("Got BAM header for iteration {}".format(iteration))
    return stdout[0]


def get_locus_names_from_bam_header(log, header, iteration):
    hs = header[0].split("\n")
    locus_names = []
    for item in hs:
        if item.startswith("@SQ"):
            locus_names.append(item.split("\t")[1].lstrip("SN:"))
    locus_names.sort()
    return locus_names


def samtools_get_locus_names_from_bam(log, bam, iteration):
    #pdb.set_trace()
    cmd1 = [
        get_user_path("executables", "samtools"),
        "view",
        bam
    ]
    cmd2 = [
        get_user_path("executables", "gawk"),
        '{print $3}',
    ]
    proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=subprocess.PIPE)
    proc1.stdout.close()
    stdout = proc2.communicate()
    # return unique list of locus names
    locus_names = list(set(stdout[0].split("\n")))
    locus_names.sort()
    # make sure empty is removed
    locus_names.remove('')
    log.info("Recovered {} loci for iteration {}".format(len(locus_names), iteration))
    return locus_names


def samtools_split_sam(sample, iteration, sample_dir_iter, sample_dir_iter_locus, locus, clean, only_single_locus):
    #pdb.set_trace()
    sam_out_fname = os.path.join(sample_dir_iter, "loci", "iter-{}.reduce.{}.bam".format(iteration, locus))
    #sam_out_fname = os.path.join(sample_dir_iter_locus, '{}.sam'.format(locus))
    # make the output dir for the locus specific files
    os.makedirs(sample_dir_iter_locus)
    # split the reduced files into properly paired and singleton reads
    bam_out_fname_paired = os.path.join(sample_dir_iter_locus, '{}.paired.bam'.format(locus))
    # -f 2 -F 2048 gets properly paired, non-supplementary alignments
    cmd2 = [
        get_user_path("executables", "samtools"),
        "view",
        "-f",
        "2",
        "-F",
        "2048",
        "-b",
        sam_out_fname,
        "-o",
        bam_out_fname_paired
    ]
    proc2 = subprocess.Popen(cmd2)
    stdout = proc2.communicate()
    # sort the paired bam
    bam_out_fname_paired_sorted = os.path.join(sample_dir_iter_locus, '{}.paired.sorted.bam'.format(locus))
    cmd1 = [
        get_user_path("executables", "samtools"),
        "sort",
        "-n",
        bam_out_fname_paired,
        "-o",
        bam_out_fname_paired_sorted
    ]
    proc1 = subprocess.Popen(cmd1)
    stdout = proc1.communicate()
    bam_out_fname_singleton = os.path.join(sample_dir_iter_locus, '{}.singleton.bam'.format(locus))
    cmd3 = [
        get_user_path("executables", "samtools"),
        "view",
        "-f",
        "8",
        "-b",
        sam_out_fname,
        "-o",
        bam_out_fname_singleton
    ]
    proc3 = subprocess.Popen(cmd3)
    stdout = proc3.communicate()
    if clean:
        os.remove(sam_out_fname)
        os.remove(bam_out_fname_paired)
    return bam_out_fname_paired_sorted, bam_out_fname_singleton


def faster_split_bam(log, sorted_reduced_bam, sample_dir_iter, iteration):
    start_time = time.time()
    sample_dir_iter_locus_temp = os.path.join(sample_dir_iter, "loci", "temp")
    # make a temp dir in locus folder in which to store locus-specific SAM data
    os.makedirs(sample_dir_iter_locus_temp)
    os.chdir(sample_dir_iter_locus_temp)
    cmd1 = [
        get_user_path("executables", "samtools"),
        "view",
        sorted_reduced_bam  
    ]
    cmd2 = [
        get_user_path("executables", "grep"),
        "-v",
        "^@"
    ]
    cmd3 = [
        get_user_path("executables", "gawk"),
        "-F\t",
        '{print > $3}'
    ]
    proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=subprocess.PIPE)
    proc3 = subprocess.Popen(cmd3, stdin=proc2.stdout, stdout=subprocess.PIPE)
    proc1.stdout.close()
    proc2.stdout.close()
    stdout = proc3.communicate()
    if proc3.returncode is not 0:
        raise IOError("Splitting BAM file has failed")
    else:
        os.chdir(sample_dir_iter)
    end_time = time.time()
    time_delta_sec = round(end_time - start_time, 3)
    log.info("\tSplit SAMs took {} seconds".format(time_delta_sec))
    return sample_dir_iter_locus_temp


def reheader_split_sams(log, sample_dir_iter, sample_dir_iter_locus_temp, header, locus_names):
    start_time = time.time()
    for locus in locus_names:
        sample_dir_iter_locus = os.path.join(sample_dir_iter, "loci", locus)
        os.makedirs(sample_dir_iter_locus)
        with open(os.path.join(sample_dir_iter_locus, "{}.sam".format(locus)), 'w') as outfile:
            with open(os.path.join(sample_dir_iter_locus_temp, locus), 'r') as temp_sam:
                outfile.write(header)
                outfile.write(temp_sam.read())
    end_time = time.time()
    time_delta_sec = round(end_time - start_time, 3)
    log.info("\tReheadering SAMs took {} seconds".format(time_delta_sec))