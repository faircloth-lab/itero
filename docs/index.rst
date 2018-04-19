.. itero documentation master file, created by
   sphinx-quickstart on Thu Apr 19 14:29:18 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. include:: global.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

itero: guided contig assembly for target enrichment data
========================================================

Release v\ |version|. (:ref:`Changelog`)

:Author: Brant C. Faircloth
:Date: |date|
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

itero_ is a "guided-assembly" workflow for target enrichment data that uses bwa_, samtools_, bedtools_, and spades_ to provided high-quality assemblies of raw reads from Illumina_ instruments.

Contributions
--------------

itero_ is open-source (see :ref:`License`) and we welcome contributions from
anyone who is interested.  Please make a pull request on github_.  The issue
tracker for itero_ is also `available on github <https://github.com/faircloth-
lab/phyluce/issues>`_.

Issues
------

If you have an issue, please ensure that you are experiencing this issue on a
supported OS (see :ref:`Installation`) using the conda_ installation of
itero_.  If possible, please submit a test case demonstrating the issue and
indicate which platform, git checkout, and phyluce version you are using.

Guide
=====

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   purpose
   installation
   running

Project info
============
.. toctree::
   :maxdepth: 1

   license
   changelog
   funding
