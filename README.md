# itero
A pipeline for iterative, guided contig assembly that integrates spades, bwa, and samtools to produce assembled contigs.  MPI and multiprocessing aware.

## Requirements
If installed using `conda`, all of the `itero` dependencies should also be installed.  These include:

### Python Modules

* numpy
* biopython
* mpi4py
* schwimmbad
* six

### 3rd-party programs

* bedtools
* bwa
* gawk
* samtools
* spades

## Installation

Easiest using `conda` and `bioconda` (**not yet implemented**).  Also installs by cloning/downloading this repository and running:

`python setup.py install`

This command should install `itero` and its `python` dependencies, but you are on your own to install the 3rd-party programs that are also needed.  The `conda` installation automates this process.

## Configuration

If you install `itero` without using `conda`, `itero` needs to know the paths to the 3rd-party binaries listed above.  You can directly specify the paths (or override the defaults), by creating an `~/.itero.conf` file with a format identical to:

```
[executables]
bedtools:$CONDA_PREFIX/bin/bedtools
bwa:$CONDA_PREFIX/bin/bwa
gawk:$CONDA_PREFIX/bin/gawk
samtools:$CONDA_PREFIX/bin/samtools
spades:$CONDA_PREFIX/bin/spades.py
```

## Execution

`itero` has both a *local* mode and an *MPI* mode.  The *local* mode is for execution on a single node, while the *MPI* mode executes individual locus assemblies in parallel using an MPI-enabled HPC cluster.  To run the program, you must first create a configuration file denoting the samples you wish you assemble. That file has the following format:

```
[reference]
/path/to/the/locus/seeds.fasta

[individuals]
taxon-one:/path/to/fastq/R1/and/R2/files/for/taxon/1/
taxon-two:/path/to/fastq/R1/and/R2/files/for/taxon/2/
taxon-three:/path/to/fastq/R1/and/R2/files/for/taxon/3/
```

You then run the *local* version using a command similar to:

```
itero assemble local --config ndna-test.conf 
	--output local
	--local-cores 16
	--iterations 6
```

This will run `itero` on a single node and will first use 16 cores to perform `bwa` alignments.  The code will then distribute locus-specific assemblies across all cores on the node (1 assembly per core; 16 in parallel).

You run the *MPI* version using a command similar to:

```
mpirun -hostfile hostfile -n 96 itero assemble mpi --config ndna-test.conf \
	--output mpi \
	--local-cores 16 \
	--iterations 6
```

If each of your nodes has 16 cores, this will first use 16 cores for the needed `bwa` alignments of reads to seeds.  The code will then distribute locus-specific assemblies across all 96 cores in your cluster (1 assembly per core; 96 in parallel).