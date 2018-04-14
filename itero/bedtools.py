#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2018 Brant Faircloth || http://faircloth-lab.org/

All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 14 April 2018 16:11 CDT (-0500)
"""

import os
import subprocess

from itero.pth import get_user_path


def bedtools_to_fastq(sample, sample_dir, bam_paired, bam_singleton, locus, clean):
    fastq_out_fname_r1 = os.path.join(sample_dir, '{}.read1.fastq'.format(locus))
    fastq_out_fname_r2 = os.path.join(sample_dir, '{}.read2.fastq'.format(locus))
    fastq_out_fname_s = os.path.join(sample_dir, '{}.singleton.fastq'.format(locus))
    cmd0 = [
        get_user_path("executables", "bedtools"),
        "bamtofastq",
        "-i",
        bam_paired,
        "-fq",
        fastq_out_fname_r1,
        "-fq2",
        fastq_out_fname_r2
    ]
    proc0 = subprocess.Popen(cmd0, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # stderr may contain entries when chimeric reads are present.  these are not
    # included in the output.
    stdout, stderr = proc0.communicate()
    cmd1 = [
        get_user_path("executables", "bedtools"),
        "bamtofastq",
        "-i",
        bam_singleton,
        "-fq",
        fastq_out_fname_s
    ]
    proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc1.communicate()
    fastqs = {
        1: fastq_out_fname_r1,
        2: fastq_out_fname_r2,
        's': fastq_out_fname_s
    }
    if clean:
        os.remove(bam_paired)
        os.remove(bam_singleton)
    return fastqs