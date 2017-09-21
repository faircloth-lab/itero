#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2017 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under an MIT License. Please see
LICENSE.txt for more information.

Created on 20 September 2017 17:22 CDT (-0500)


Parts of the code below is derived from embarrassingly_parallel.py
(https://gist.github.com/krischer/2c7b95beed642248487a) by Lion Krischer which
is available under the MIT License.


The MIT License (MIT)

Copyright (c) 2015 Lion Krischer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""


import os
import sys
import numpy

from itero import samtools_split_bam, bedtools_to_fastq, spades_paired_end_assembly, initial_assembly
from mpi4py import MPI


def split(container, count):
    """
    Simple function splitting a container into equal length chunks.
    Order is not preserved but this is potentially an advantage depending on
    the use case.
    """
    return [container[_i::count] for _i in range(count)]


COMM = MPI.COMM_WORLD
if COMM.rank == 0:
    # take parameters coming in as arguments - this script is called from the
    # main loop because it makes things run much cleaner (also gives us some
    # options relative to running both mutltprocessing and mpi)
    iteration = sys.argv[1]
    sample = sys.argv[2]
    sample_dir_iter = sys.argv[3]
    sorted_reduced_bam = sys.argv[4]
    if sys.argv[5] == 'True':
        clean = True
    else:
        clean = False
    if sys.argv[6] == 'True':
        multiple_hits = True
    else:
        multiple_hits = False
    locus_file = os.path.join(sample_dir_iter, sys.argv[6])
    locus_names = numpy.genfromtxt(locus_file, delimiter=',', dtype=str)
    # generate a work tuple as in the multiprocessing version of the code
    work = [(iteration, sample, sample_dir_iter, sorted_reduced_bam, locus_name, clean, multiple_hits) for locus_name in locus_names]
    jobs = split(work, COMM.size)
else:
    jobs = None

# Scatter jobs across cores.
jobs = COMM.scatter(jobs, root=0)

# Now each rank just does its jobs. We don't need a results lists because we get
# that in the main loop - basically we just need the dumbest parallelization
# possible. Make sure to not use super big objects in there as they will be
# pickled to be exchanged over MPI.

for job in jobs:
    initial_assembly(job)
