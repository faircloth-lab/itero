#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2018 Brant Faircloth || http://faircloth-lab.org/

All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 14 April 2018 16:09 CDT (-0500)
"""

from itero.helpers import FullPaths, CreateDir, is_dir, is_file

def get_common_arguments(parser):
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
    return parser