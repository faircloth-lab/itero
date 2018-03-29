
import os
import subprocess

from itero.pth import get_user_path


def samtools_index(log, sample, sample_dir, bam, iteration=0):
    log.info("Indexing BAM for {}".format(sample))
    cmd = [
        get_user_path("executables", "samtools"),
        "index",
        bam
    ]
    samtools_out_fname = os.path.join(sample_dir, 'iter-{}.samtools-idx.log'.format(sample))
    with open(samtools_out_fname, 'w') as samtools_out:
        proc = subprocess.Popen(cmd, stdout=samtools_out, stderr=subprocess.STDOUT)
        proc.communicate()


def samtools_reduce(log, sample, sample_dir, bam, iteration=0):
    #pdb.set_trace()
    log.info("Reducing BAM for {}, iteration {}".format(sample, iteration))
    bam_out_fname = os.path.join(sample_dir, 'iter-{}.reduce.bam'.format(iteration))
    cmd = [
        get_user_path("executables", "samtools"),
        "view",
        "-F",
        "4",
        "-bq",
        "1",
        bam,
        "-o",
        bam_out_fname
    ]
    samtools_out_fname = os.path.join(sample_dir, 'iter-{}.reduce.log'.format(iteration))
    with open(samtools_out_fname, 'w') as samtools_out:
        proc = subprocess.Popen(cmd, stdout=samtools_out, stderr=subprocess.STDOUT)
        proc.communicate()
    return bam_out_fname


def samtools_sort(log, sample, sample_dir, bam, iteration=0):
    #pdb.set_trace()
    bam_out_fname = os.path.join(sample_dir, 'iter-{}.reduce.sorted.bam'.format(iteration))
    cmd1 = [
        get_user_path("executables", "samtools"),
        "sort",
        bam,
        "-o",
        bam_out_fname
    ]
    samtools_out_fname = os.path.join(sample_dir, 'iter-{}.sort.log'.format(iteration))
    with open(samtools_out_fname, 'w') as samtools_out:
        proc = subprocess.Popen(cmd1, stdout=samtools_out, stderr=subprocess.STDOUT)
        proc.communicate()
    return bam_out_fname


def samtools_get_locus_names_from_bam(log, bam, iteration):
    #pdb.set_trace()
    cmd1 = [
        get_user_path("executables", "samtools"),
        "view",
        bam
    ]
    cmd2 = [
        get_user_path("executables", "gawk"),
        '{print $3}',
    ]
    proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=subprocess.PIPE)
    proc1.stdout.close()
    stdout = proc2.communicate()
    # return unique list of locus names
    locus_names = list(set(stdout[0].split("\n")))
    locus_names.sort()
    # make sure empty is removed
    locus_names.remove('')
    log.info("Recovered {} loci for iteration {}".format(len(locus_names), iteration))
    return locus_names


def samtools_split_bam(sample, sample_dir, bam, locus, clean, only_single_locus):
    bam_out_fname = os.path.join(sample_dir, '{}.bam'.format(locus))
    if only_single_locus:
        bam_out_fname = bam
    else:
        cmd0 = [
            get_user_path("executables", "samtools"),
            "view",
            "-b",
            bam,
            locus,
            "-o",
            bam_out_fname
        ]
        proc0 = subprocess.Popen(cmd0)
        stdout = proc0.communicate()
    # split the reduced files into properly paired and singleton reads
    bam_out_fname_paired = os.path.join(sample_dir, '{}.paired.bam'.format(locus))
    # -f 2 -F 2048 gets properly paired, non-supplementary alignments
    cmd2 = [
        get_user_path("executables", "samtools"),
        "view",
        "-f",
        "2",
        "-F",
        "2048",
        "-b",
        bam_out_fname,
        "-o",
        bam_out_fname_paired
    ]
    proc2 = subprocess.Popen(cmd2)
    stdout = proc2.communicate()
    # sort the paired bam
    bam_out_fname_paired_sorted = os.path.join(sample_dir, '{}.paired.sorted.bam'.format(locus))
    cmd1 = [
        get_user_path("executables", "samtools"),
        "sort",
        "-n",
        bam_out_fname_paired,
        "-o",
        bam_out_fname_paired_sorted
    ]
    proc1 = subprocess.Popen(cmd1)
    stdout = proc1.communicate()
    bam_out_fname_singleton = os.path.join(sample_dir, '{}.singleton.bam'.format(locus))
    cmd3 = [
        get_user_path("executables", "samtools"),
        "view",
        "-f",
        "8",
        "-b",
        bam_out_fname,
        "-o",
        bam_out_fname_singleton
    ]
    proc3 = subprocess.Popen(cmd3)
    stdout = proc3.communicate()
    if clean:
        os.remove(bam_out_fname)
        os.remove(bam_out_fname_paired)
    return bam_out_fname_paired_sorted, bam_out_fname_singleton