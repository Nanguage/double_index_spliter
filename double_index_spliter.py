import re
from collections import namedtuple
import logging
import os
from os.path import join, exists
import gzip
import sys
from sys import stderr
from io import TextIOWrapper

import click
from cutadapt.align import Aligner
from Bio import SeqIO


log = logging.getLogger(__name__)


class Index(object):
    """ Denote a kind of index(a,b) combination """
    def __init__(self, name, index_a, index_b, phred=33):
        self.name = name
        self.index_a = index_a
        self.index_b = index_b
        self.phred=phred

    def fopen(self, path="./", suffix=".fq"):
        self.file = open_(join(path, self.name+suffix), 'w')
        self.writer = fastq_writer(self.file, self.phred)
        self.writer.write_header()

    def __del__(self):
        if hasattr(self, 'file'):
            try:
                self.writer.write_footer()
            except AssertionError as e:
                msg = "Index combination: \"%s\" %s not found any reads" % \
                      (self.name, str(self))
                log.warning(msg)
            self.file.close()
    
    def __str__(self):
        if self.index_a and self.index_b:
            return self.index_a + "+" + self.index_b
        else:
            return "<unmatch>"
    
    def _aligner(self, mismatch):
        if not hasattr(self, "_aligner_a"):
            max_err_r = (float(mismatch[0])/len(self.index_a),
                         float(mismatch[1])/len(self.index_b))
            self._aligner_a = Aligner(self.index_a, max_err_r[0])
            self._aligner_b = Aligner(self.index_b, max_err_r[1])
            self._aligner_a.min_overlap = len(self.index_a) - mismatch[0]
            self._aligner_b.min_overlap = len(self.index_b) - mismatch[1]
        return self._aligner_a, self._aligner_b

    def match(self, idx_a, idx_b, mismatch):
        alr_a, alr_b = self._aligner(mismatch)
        if len(idx_a) != len(self.index_a) or \
           len(idx_b) != len(self.index_b):
           return False
        res_a = alr_a.locate(idx_a)
        if not res_a:
            return False
        res_b = alr_b.locate(idx_b)
        if not res_b:
            return False
        return True


def read_indexes(f_path, phred):
    """
    :f_path: path to indexes config file
    """
    indexes = []
    with open(f_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            items = line.split()
            n, a, b = items[:3]
            indexes.append(Index(name=n, index_a=a, index_b=b, phred=phred))
    return indexes


def get_mismatch(mismatch_str):
    # get mismatch on index-a and index-b
    # by parse mismatch_str
    mismatch = mismatch_str
    if "," in mismatch:
        mis_ = mismatch.split(",")
        mis_a, mis_b = int(mis_[0]), int(mis_[1])
    else:
        mis_a = mis_b = int(mismatch)
    return mis_a, mis_b


def open_(fname, mode):
    """ open file according to it's suffix """
    assert 'b' not in mode
    if fname.endswith(".gz"):
        if sys.version_info > (3, 0, 0):
            return TextIOWrapper(gzip.open(fname, mode+'b'))
        else:
            return gzip.open(fname, mode+'b')
    else:
        return open(fname, mode)


def fastq_iter(file_in, phred):
    """ return a fastq iterator """
    if str(phred) == '33':
        fastq_iter = SeqIO.parse(file_in, 'fastq')
    else:
        fastq_iter = SeqIO.parse(file_in, 'fastq-illumina')
    return fastq_iter


def fastq_writer(file_out, phred):
    """ return a fastq writer """
    if str(phred) == '33':
        writer = SeqIO.QualityIO.FastqPhredWriter(file_out)
    else:
        writer = SeqIO.QualityIO.FastqIlluminaWriter(file_out)
    return writer


def extract_index(id_str, a_b_spliter='+'):
    """ extract index-a and index-b from fastq id string """
    idx_str = id_str.split(":")[-1]
    idx_a, idx_b = idx_str.split(a_b_spliter)
    return idx_a, idx_b


def process_all(input_fq, all_indexes, unmatched, mismatch, phred):
    with open_(input_fq, 'r') as f:
        fq_iter = fastq_iter(f, phred)
        for rec in fq_iter:
            idx_a, idx_b = extract_index(rec.description)
            for idx in all_indexes:
                if idx.match(idx_a, idx_b, mismatch):
                    idx.writer.write_record(rec)
                    break
            else:
                unmatched.writer.write_record(rec)


@click.command(name="double_index_spliter")
@click.argument("input_fastq")
@click.argument("index_file")
@click.option("--mismatch", "-m",
    default="2,2",
    help="mismatch threshold of index-a and index-b, i.e. -m 2,3 means "
         "2 bases mismatch on index-a, and 3 bases mismatch on b")
@click.option("--outdir", "-O",
    default="./",
    help="path to output splited fastq files.")
@click.option("--gzip/--no-gzip", "-z",
    default=False,
    help="compress output fastq files with gzip.")
@click.option("--phred",
    default=33,
    help="encode of fastq quality string.")
def main_(input_fastq, index_file, mismatch, outdir, gzip, phred):
    mis_a, mis_b = mismatch = get_mismatch(mismatch)
    log.info("mismatch threshold:\nindex-a: {}\nindex-b: {}".format(mis_a, mis_b))

    log.info("parsing index config file: %s"%index_file)
    all_idx = read_indexes(index_file, phred)
    unmatch = Index("unmatched", index_a=None, index_b=None, phred=phred)

    if not exists(outdir):
        os.mkdir(outdir)
    # touch all sub fq file
    suffix = ".fq" if not gzip else ".fq.gz"
    for idx in all_idx:
        idx.fopen(path=outdir, suffix=suffix)
        msg = "splited index: {} 's fastq will store to {}".format(
            str(idx), join(outdir, idx.name)+suffix )
        log.info(msg)
    # touch unmatched fq file
    unmatch.fopen(path=outdir, suffix=suffix)
    
    process_all(input_fastq, all_idx, unmatch, mismatch, phred)


if __name__ == "__main__":
    log.setLevel(logging.DEBUG)
    hdr = logging.StreamHandler(stream=stderr)
    fmt = logging.Formatter(
        fmt="[%(levelname)s][%(asctime)s] %(message)s",
        datefmt="%m/%d/%y %H:%M:%S")
    hdr.setFormatter(fmt)
    log.addHandler(hdr)
    main_()
