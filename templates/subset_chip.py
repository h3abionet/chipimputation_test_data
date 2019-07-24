#!/usr/bin/env python
"""
Reads chip file, chunk size
Returns file with chromosome chunk_start chunk_end
"""

import argparse,sys,time
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--chip", default="${chip}", help="")
parser.add_argument("--subset_map", default="${subset_map}", help="")
parser.add_argument("--out_chip", default="${out_chip}", help="")
args = parser.parse_args()

def subset_map(chip, subset_map, out_chip):
    '''
    Return: max and min base for each chromosome
    '''
    data = {}
    out = open(out_chip, 'w')

    for line in open(chip):
        dat = line.strip().split(',')
        try:
            chr = dat[9]
            pos = dat[10]
            chr_pos = chr+"_"+pos
            data[chr_pos] = '\t'.join([chr, pos, pos])+'\\n'
        except:
            continue
    for line in open(subset_map):
        dat = line.strip().split()
        try:
            chr = dat[0]
            pos = dat[1]
            chr_pos = chr+"_"+pos
            if chr_pos in data:
                out.writelines(data[chr_pos])
        except:
            continue

    out.close()


if __name__ == '__main__':
    subset_map(args.chip, args.subset_map, args.out_chip)

