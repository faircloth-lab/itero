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
import argparse

import itero
from itero.cli import main_help
from itero.cli import main_assemble
from itero.cli import main_check


def main():
    # print same output as help if no arguments
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    # setup main program args
    p = argparse.ArgumentParser(
        description="itero is a software package for iterative, guided assembly of "
                    "target enrichment data."
    )
    p.add_argument(
        "-V", "--version",
        action="version",
        version="itero {}".format(itero.__version__)
    )

    sub_parsers = p.add_subparsers(
        metavar="command",
        dest="cmd",
    )

    main_help.configure_parser(sub_parsers)
    main_assemble.configure_parser(sub_parsers)
    main_check.configure_parser(sub_parsers)
    
    try:
        import argcomplete
        argcomplete.autocomplete(p)
    except ImportError:
        pass
    except AttributeError:
        # On Python 3.3, argcomplete can be an empty namespace package when
        # argcomplete is not installed. Not sure why, but this fixes it.
        pass

    args = p.parse_args()

    try:
        args.func(args, p)
    except RuntimeError as e:
        sys.exit("Error: %s" % e)
    except Exception as e:
        raise

    sys.exit(0)


if __name__ == "__main__":
    main()