#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2018 Brant Faircloth || http://faircloth-lab.org/

All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 21 December 2018 10:45 CST (-0600)
"""

import os
import time
import subprocess
import shutil

from itero.pth import get_user_path

import pdb

def split_bam(log, sorted_reduced_bam, sample_dir_iter, iteration):
    start_time = time.time()
    sample_dir_iter_locus = os.path.join(sample_dir_iter, "loci")
    # make a temp dir in locus folder in which to store locus-specific SAM data
    os.makedirs(sample_dir_iter_locus)
    os.chdir(sample_dir_iter_locus)
    # copy the sorted BAM into this directory
    temp_bam_pth = os.path.join(sample_dir_iter_locus, os.path.basename(sorted_reduced_bam))
    shutil.copyfile(sorted_reduced_bam, temp_bam_pth)
    #pdb.set_trace()
    cmd1 = [
        get_user_path("executables", "bamtools"),
        "split",
        "-in",
        temp_bam_pth,
        "-reference",
        "-refPrefix",
        ""
    ]
    proc1 = subprocess.Popen(cmd1)
    stdout = proc1.communicate()
    if proc1.returncode is not 0:
        raise IOError("Splitting BAM file has failed")
    else:
        os.chdir(sample_dir_iter)
    os.remove(temp_bam_pth)
    end_time = time.time()
    time_delta_sec = round(end_time - start_time, 3)
    log.info("\tSplit BAMs took {} seconds".format(time_delta_sec))
    return sample_dir_iter_locus