#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2018 Brant Faircloth || http://faircloth-lab.org/

All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 14 April 2018 16:12 CDT (-0500)
"""

import os
import sys
import ConfigParser

#import pdb


def get_user_path(program, binary, package_only=False):
    config = ConfigParser.ConfigParser()
    # make case sensitive
    config.optionxform = str
    if package_only:
        config.read(os.path.join(sys.prefix, 'config/itero.conf'))
    else:
        config.read([os.path.join(sys.prefix, 'config/itero.conf'), os.path.expanduser('~/.itero.conf')])
    # ensure program is in list
    pth = config.get(program, binary)
    # expand path as necessary - replace CONDA variable placeholder
    # with sys.prefix, otherwise default to normal path expansion
    if pth.startswith("$CONDA_PREFIX"):
        expand_pth = pth.replace("$CONDA_PREFIX", sys.prefix)
    else:
        expand_pth = os.path.abspath(os.path.expanduser(os.path.expandvars(pth)))
    return expand_pth
