#!/usr/bin/env python
"""
"""

import argparse,sys,time
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--fasta_in", default="${fasta_in}", help="")
parser.add_argument("--fasta_out", default="${fasta_out}", help="")
args = parser.parse_args()

def fix_fasta_header(fasta_in, fasta_out):
    '''
    '''
    nline = 1
    out = open(fasta_out, 'w')
    for line in open(fasta_in):
        if line.startswith(">"):
            line = line[1:].strip().split(":")
            try:
                chr = line[0]
                pos = line[1].split("-")
                start = pos[1]
                end = str(int(pos[1])+1)
                line = ">{3} dna:chromosome chromosome:GRCh37:{0}:{1}-{2}:1".format(chr, start, end, nline)+"\n"
                nline += 1
            except:
                continue
        out.writelines(line)
    out.close()


if __name__ == '__main__':
    fix_fasta_header(args.fasta_in, args.fasta_out)

