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
import subprocess

from itero.pth import get_user_path

#import pdb


def spades_version():
    cmd = [get_user_path("executables", "spades"), '--version']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = proc.communicate()
    return stdout.split("\n")[0].split(' ')[1]


def spades_paired_end_assembly(iteration, sample, sample_dir, fastqs, locus, clean):
    assembly_out_fname = os.path.join(sample_dir, '{}-assembly'.format(locus))
    # go ahead and assemble without error correction, for speed.
    # explcitly set threads = 1
    cmd1 = [
        get_user_path("executables", "spades"),
        "-t",
        "1",
        "-1",
        fastqs[1],
        "-2",
        fastqs[2],
        "-s",
        fastqs['s'],
        "-k",
        "33",
        "--cov-cutoff",
        "5",
        "-o",
        assembly_out_fname
    ]
    # turn off error correction for non-final rounds, turn on error-correction
    # for final round and also use --careful assembly option (both of these are
    # slower)
    if not iteration == 'final':
        cmd1.append("--only-assembler")
    if iteration == 'final':
        cmd1.append("--careful")
    # spades creates its own log file in the assembly dir - redirect to /dev/null
    fnull_file = open(os.devnull, 'w')
    proc = subprocess.Popen(cmd1, stdout=fnull_file, stderr=subprocess.STDOUT)
    stdout, stderr = proc.communicate()
    return assembly_out_fname