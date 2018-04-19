.. include:: global.rst

.. _Installation:

************
Installation
************

itero_ uses a number of Python_ tools that allow it to assemble raw reads into contigs.  itero_ also wraps a number of third-party programs.  These include:

Python Modules
===============

* numpy_
* biopython_
* mpi4py_
* schwimmbad_
* six_

3rd-party programs
===================

* bedtools_
* bwa_
* gawk_
* samtools_
* spades_

To ensure that these dependencies are easy to install, we have created a conda_ package for itero_ that is available as part of bioconda_.  This is the easiest way to get itero up and running on your system.  itero_ can also be run outside of conda_, and we include some installation suggestions for these types of systems, below.  However, because many HPC systems are configured different, we cannot provide extended support for itero_ on HPC platforms.

.. note:: We build and test the binaries available through conda_ using
    64-bit operating systems that include:

    - Apple OSX 10.9.x
    - CentOS 7.x

Why conda?
==========

It may seem odd to impose a particular disitribution on users, and we largely
agree.  However, conda_ makes it very easy for us to distribute both Python_ and
non-Python packages, setup identical environments
across very heterogenous platforms (linux, osx), make sure all the `$PATHs` are
correct, and have things run largely as expected. Using conda_ has several other
benefits, including environment separation similar to virtualenv_.

In short, using conda_ gets us as close to a "one-click" install that we will probably
ever get.

Install Process (using conda/bioconda)
======================================

.. attention:: We do not support itero_ on Windows.

.. note:: We build and test the binaries available through using
   64-bit operating systems that include:

   - Apple OSX 10.9.x
   - CentOS 7.x

The installation process is a 3-step process.  You need to:

#. Install conda_ (either anaconda_ or miniconda_)
#. Configure conda_ to use bioconda_
#. Install itero_

Installing itero_ using conda_ will install all of the required binaries, libraries, and
Python_ dependencies.

Install Anaconda or miniconda
-----------------------------

You first need to install anaconda_ or miniconda_.  Which
one you choose is up to you, your needs, how much disk space you have, and if
you are on a fast/slow connection.

.. attention:: You can easily install anaconda_ or miniconda_ in your $HOME,
    although you should be aware that this setup can cause problems in some
    HPC setups.

.. tip:: Do I want anaconda_ or miniconda_?
    :class: admonition tip

    The major difference between the two python distributions is that anaconda_
    comes with many, many packages pre-installed, while miniconda_ comes with
    almost zero packages pre-installed.  As such, the beginning anaconda_
    distribution is roughly 500 MB in size while the beginning miniconda_
    distribution is 15-30 MB in size.


anaconda
^^^^^^^^

Follow the instructions here for your platform:
http://docs.continuum.io/anaconda/install.html

miniconda
^^^^^^^^^

Find the correct `miniconda-x.x.x` file for your platform from
http://repo.continuum.io/miniconda/ and download that file.  Be sure you **do
not** get one of the packages that has a name starting with `miniconda3-`. When
that has completed, run one of the following::

    bash Miniconda-x.x.x-Linux-x86_64.sh  [linux]
    bash Miniconda-x.x.x-MacOSX-x86_64.sh [osx]

.. note:: Once you have installed Miniconda, we will refer to it as **anaconda**
   throughout the remainder of this documentation.

Checking your `$PATH`
^^^^^^^^^^^^^^^^^^^^^

Regardless of whether you install anaconda_ or miniconda_, you need to check
that you've installed the package correctly.  To ensure that the correct
location for anaconda_ or miniconda_ are added to your $PATH (this occurs
automatically on the $BASH shell), run the following::

    $ python -V

The output should look similar to (`x` will be replaced by a version)::

    Python 2.7.x :: Anaconda x.x.x (x86_64)

Notice that the output shows we're using the `Anaconda x.x.x` version of
Python_. If you do not see the expected output (or something similar), then you
likely need to edit your $PATH variable to add anaconda_ or miniconda_.

The easiest way to edit your path, if needed is to open ``~/.bashrc`` with a
text editor (if you are using ZSH, this will be ``~/.zshrc``) and add, as the
last line::

    export PATH=$HOME/path/to/conda/bin:$PATH

where ``$HOME/path/to/conda/bin`` is the location of anaconda/miniconda on your
system (usually ``$HOME/anaconda/bin`` or ``$HOME/miniconda/bin``).

.. warning:: If you have previously set your ``$PYTHONPATH`` elsewhere in your
   configuration, it may cause problems with your anaconda_ or miniconda_
   installation of phyluce_.  The solution is to remove the offending library
   (-ies) from your ``$PYTHONPATH``.

Configure Bioconda
------------------

Once you have installed anaconda_ (or miniconda_), you need to configure conda_ to use the bioconda_ channel.  More information on this process can be found on the bioconda_ website, but the gist of the process is that you need to run::

    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

Install itero_
--------------

Once bioconda_ is installed, you should be able to install itero_ by running::
    
    conda install itero

This should install everything that you need to run the program.

Test itero_ install
-------------------

You can check to make sure all of the binaries are installed correctly by running::

    itero check binaries


Install Process (Alternative / HPC)
===================================

On some systems (particularly HPC systems), conda_ can cause problems.  You can install the itero_ package the "old" way by downloading the package tarball (https://github.com/faircloth-lab/itero/releases) and running::

    python setup.py install

in the main directory.  This should install all of the Python_ dependencies, **but you still need to install and configure the 3rd-party dependencies**.

.. attention:: You will need to install 3rd-party dependencies on your own
  if you are using the `python setup.py install` method of installing itero_

You can build and install these dependencies where you like.  To configure itero_ to use the dependencies you have build and installed, you need to create a ``$HOME/.itero.conf`` that gives the paths to each program and looks like:

.. code-block:: text

    [executables]
    bedtools:/path/to/bin/bedtools
    bwa:/path/to/bin/bwa
    gawk:/path/to/bin/gawk
    samtools:/path/to/bin/samtools
    spades:/path/to/bin/spades.py


Test itero_ install
-------------------

You can check to make sure all of the binaries are installed correctly by running::

    itero check binaries