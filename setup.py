#!/usr/bin/env python
# encoding: utf-8

import os
from setuptools import setup
from itero import __version__


def package_files(directory):
    """package up files.  from stackoverflow.com
    how-to-add-package-data-recursively-in-python-setup-py"""
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths


setup(
    name='itero',
    version=__version__,
    description='a strongly opinionated, iterative, guided assembly pipeline',
    url='https://github.com/faircloth-lab/itero',
    author='Brant C. Faircloth',
    author_email='borg@faircloth-lab.org',
    license='MIT',
    platforms='any',
    packages=[
        'itero', 'itero.assemble', 'itero.check', 'itero.cli'
    ],
    data_files=[('config', ['itero/config/itero.conf'])],
    scripts=['bin/itero'],
    package_data={'': package_files('itero/tests')},
    install_requires=['numpy', 'schwimmbad', 'biopython', 'six', 'mpi4py'],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)