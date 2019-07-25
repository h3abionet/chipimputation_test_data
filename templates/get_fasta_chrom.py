#!/usr/bin/env python
"""
"""

import argparse,sys,time

parser = argparse.ArgumentParser()
parser.add_argument("--fasta_in", default="${fasta_in}", help="")
parser.add_argument("--chrm", default="${chrm}", help="")
parser.add_argument("--fasta_out", default="${fasta_out}", help="")
args = parser.parse_args()

def get_fasta_chrom(fasta_in, chrm, fasta_out):
    """
    :param fasta_in:
    :param chrm:
    :param fasta_out:
    :return:
    """
    out = open(fasta_out, 'w')
    for line in open(fasta_in):
        if line.startswith(">") and 'chromosome:GRCh37:'+str(chrm) in line:
            cond = True
        else:
            if line.startswith(">"):
                cond = False
        if cond:
            out.writelines(line)
        else:
            continue

    out.close()

if __name__ == '__main__':
    get_fasta_chrom(args.fasta_in, args.chrm, args.fasta_out)

