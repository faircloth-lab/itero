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
import time
import logging


import pdb

def setup_logging(args):
    import __main__ as main
    import __init__ as init
    my_name = os.path.basename(os.path.splitext(main.__file__)[0])
    log = logging.getLogger(my_name)
    console = logging.StreamHandler(sys.stdout)
    curr_time = int(time.time())
    if args.log_path is not None:
        logfile = logging.FileHandler(os.path.join(args.log_path, "{}.log".format(my_name, curr_time)))
    else:
        logfile = logging.FileHandler("{}-{}.log".format(my_name, curr_time))
    if args.verbosity == "INFO":
        log.setLevel(logging.INFO)
        console.setLevel(logging.INFO)
        logfile.setLevel(logging.INFO)
    if args.verbosity == "WARN":
        log.setLevel(logging.WARN)
        console.setLevel(logging.WARN)
        logfile.setLevel(logging.WARN)
    if args.verbosity == "CRITICAL":
        log.setLevel(logging.CRITICAL)
        console.setLevel(logging.CRITICAL)
        logfile.setLevel(logging.CRITICAL)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logfile.setFormatter(formatter)
    log.addHandler(console)
    log.addHandler(logfile)
    text = " Starting {} ".format(my_name)
    log.info(text.center(65, "="))
    log.info("Version: {}".format(init.__version__))
    for arg, value in sorted(vars(args).items()):
        log.info("Argument --{}: {}".format(arg, value))
    return log, my_name
