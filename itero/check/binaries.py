#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2018 Brant Faircloth || http://faircloth-lab.org/

All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 18 April 2018 15:35 CDT (-0500)
"""

from collections import OrderedDict

import numpy
import Bio
import schwimmbad

from itero import bwa
from itero import samtools
from itero import bedtools
from itero import spades

# import pdb


def main(args, parser):
    versions = (
        ("bwa", bwa.bwa_version()),
        ("samtools", samtools.samtools_version()),
        ("bedtools", bedtools.bedtools_version()),
        ("spades", spades.spades_version()),
        ("schwimmbad", schwimmbad.__version__),
        ("numpy", numpy.__version__),
        ("biopython", Bio.__version__)
    )

    for items in versions:
        print "{:<12}: {}".format(items[0], items[1])