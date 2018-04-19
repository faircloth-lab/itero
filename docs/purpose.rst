.. include:: global.rst

Purpose
=======

itero_ is a program what wraps a workflow for guided- or reference-based assembly of target enrichment data.  This approach to assembly is also "iterative", meaning that the assembly proceeds through several, iterations (hence "itero").  I wrote itero_ for a variety of reasons:

* "traditional" DNA assembly programs performed poorly with target enrichment data (from UCE loci)
* existing DNA assembly approaches had relatively high assembly error
* existing guided assembly programs were hard to install and run
* some existing guided assembly programs were **slow**

itero_ attempts to fix some of these problems.  At its heart, itero_ uses an input file of "seeds", against which it will try to assemble raw-read data from Illumina_ instruments.  Alignment of reads-to-seeds uses bwa_, the BAM file is split with samtools_ and bedtools_, and locus-specific reads are then assembled using spades_ (with error correction turned on during the final round). Then, the entire process repeats itself.

To increase assembly speed, itero_ takes advantage of multiple cores (on single nodes) using python_ multiprocessing and MPI (on HPC systems) using the excellent schwimmbad_ library.

Who wrote this?
---------------

This documentation was written primarily by Brant Faircloth
(http://faircloth-lab.org). Brant is also responsible for the development of
most of the itero_ code.  Bugs within the code are usually his.


How do I report bugs?
----------------------

To report a bug, please post an issue to
https://github.com/faircloth-lab/itero/issues.  Please also ensure that you
are using one of the "supported" platforms:

- Apple OSX 10.9.x
- CentOS 7.x

and that you have installed itero_ and dependencies using conda_, as described
in the :ref:`Installation` section.