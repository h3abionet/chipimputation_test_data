#!/usr/bin/env python
"""
Reads chip file, chunk size
Returns file with chromosome chunk_start chunk_end
"""

import argparse,sys,time,random

parser = argparse.ArgumentParser()
parser.add_argument("--subset_map", default="${subset_map}", help="")
parser.add_argument("--out_chip", default="${out_chip}", help="")
args = parser.parse_args()

def subset_map(subset_map, out_chip):
    '''
    Return: max and min base for each chromosome
    '''
    data = []
    out = open(out_chip, 'w')
    for line in open(subset_map):
        dat = line.strip().split()
        try:
            chr = dat[0]
            pos = dat[1]
            data.append('\t'.join([chr, pos, pos])+'\\n')
        except:
            continue
    random.shuffle(data)
    for dat in data[:5*len(data)/100]:
        out.writelines(dat)
    out.close()


if __name__ == '__main__':
    subset_map(args.subset_map, args.out_chip)

