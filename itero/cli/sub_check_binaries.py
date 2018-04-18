#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2018 Brant Faircloth || http://faircloth-lab.org/

All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 14 April 2018 16:10 CDT (-0500)
"""

from __future__ import absolute_import

from itero.check import binaries


descr = "Check to ensure binaries are installed and configured."


def configure_parser(sub_parsers):
    sp = sub_parsers.add_parser(
        'binaries',
        description=descr,
        help=descr
    )
    # get common arguments
    sp.set_defaults(func=run_check_binaries)


def run_check_binaries(args, parser):
    binaries.main(args, parser)