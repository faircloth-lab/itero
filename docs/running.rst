.. include:: global.rst

Running itero
=============

itero_ has both a **local** mode and an **MPI** mode.  The **local** mode is for execution on a single node, while the **MPI** mode executes individual locus assemblies in parallel using an MPI-enabled HPC cluster.  To run the program, you must first create a configuration file denoting the samples you wish you assemble. That file has the following format:

.. code-block:: text

    [reference]
    /path/to/the/locus/seeds.fasta

    [individuals]
    taxon-one:/path/to/fastq/R1/and/R2/files/for/taxon/1/
    taxon-two:/path/to/fastq/R1/and/R2/files/for/taxon/2/
    taxon-three:/path/to/fastq/R1/and/R2/files/for/taxon/3/

itero_ on a single node
-----------------------

You then run the *local* version using a command similar to:

.. code-block:: bash

    itero assemble local --config ndna-test.conf 
        --output local
        --local-cores 16
        --iterations 6


This will run `itero` on a single node and will first use 16 cores to perform `bwa` alignments.  The code will then distribute locus-specific assemblies across all cores on the node (1 assembly per core; 16 in parallel).


itero_ across MPI nodes
-----------------------

You run the *MPI* version using a command similar to:

.. code-block:: text

    mpirun -hostfile hostfile -n 96 itero assemble mpi --config ndna-test.conf \
        --output mpi \
        --local-cores 16 \
        --iterations 6

If each of your nodes has 16 cores, this will first use 16 cores for the needed `bwa` alignments of reads to seeds.  The code will then distribute locus-specific assemblies across all 96 cores in your cluster (1 assembly per core; 96 in parallel).