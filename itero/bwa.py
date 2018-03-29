

import os
import subprocess

from itero.pth import get_user_path

#import pdb


def bwa_create_index_files(log, reference):
    log.info("Running bwa indexing against {}".format(reference))
    cwd = os.getcwd()
    # move into reference directory
    os.chdir(os.path.dirname(reference))
    cmd = [get_user_path("executables", "bwa"), "index", reference]
    with open('bwa-index-file.log', 'a') as outf:
        proc = subprocess.Popen(cmd, stdout=outf, stderr=subprocess.STDOUT)
        proc.communicate()
    # mvoe back to working directory
    os.chdir(cwd)


def bwa_index_seeds(seeds, log):
    #pdb.set_trace()
    log.info("Running bwa indexing against {}".format(os.path.basename(seeds)))
    cwd = os.getcwd()
    # move into reference directory
    os.chdir(os.path.dirname(seeds))
    cmd = [get_user_path("executables", "bwa"), "index", seeds]
    with open('bwa-index-file.log', 'a') as outf:
        proc = subprocess.Popen(cmd, stdout=outf, stderr=subprocess.STDOUT)
        proc.communicate()
    # mvoe back to working directory
    os.chdir(cwd)


def bwa_mem_pe_align(log, sample, sample_dir, ref, cores, r1, r2, iteration=0):
    #pdb.set_trace()
    cmd1 = [
        get_user_path("executables", "bwa"),
        "mem",
        "-t",
        str(cores),
        ref,
        r1.pth,
        r2.pth
    ]
    cmd2 = [
        get_user_path("executables", "samtools"),
        "view",
        "-bS",
        "-"
    ]
    sampe_out_fname = os.path.join(sample_dir, 'iter-{}.pe.bwa.log'.format(iteration))
    samtools_out_fname = os.path.join(sample_dir, 'iter-{}.pe.samtools.log'.format(iteration))
    bam_out_fname = os.path.join(sample_dir, 'iter-{}.bam'.format(iteration))
    log.info("Building BAM for {}, iteration {}".format(sample, iteration))
    with open(sampe_out_fname, 'w') as sampe_out:
        with open(samtools_out_fname, 'w') as samtools_out:
            with open(bam_out_fname, 'w') as bam_out:
                proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=sampe_out)
                proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=bam_out, stderr=samtools_out)
                proc1.stdout.close()
                proc2.communicate()
    return bam_out_fname