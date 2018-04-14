#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2018 Brant Faircloth || http://faircloth-lab.org/

All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 14 April 2018 16:09 CDT (-0500)
"""

from __future__ import absolute_import

import sys

from itero.cli import sub_assemble_local
from itero.cli import sub_assemble_mpi


descr = "Assemble cleaned/trimmed sequencing reads."


def configure_parser(sub_parsers):
    if len(sys.argv) == 2:
        sys.argv.append("-h")
    
    p = sub_parsers.add_parser(
        'assemble',
        description=descr,
        help=descr
    )

    sub_parsers = p.add_subparsers(
        metavar="command",
        dest="cmd",
    )

    sub_assemble_local.configure_parser(sub_parsers)
    sub_assemble_mpi.configure_parser(sub_parsers)