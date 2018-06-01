..  _Changelog:

Changelog
=========

v1.0.x (April 2017)
-----------------

* initial version with MPI and multiprocessing capability


v1.1.x (May 2017)
-----------------

* fix error in contig checking code that could cause MPI operations to hang
* refactor BAM splitting code for hopefully faster operation
* add RAM limits on spades
* add configuration parameters to iter.conf for spades
* create unique log file for each run
