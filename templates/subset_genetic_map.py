#!/usr/bin/env python
"""
Reads map file, chunk size
Returns file with chromosome chunk_start chunk_end
"""

import argparse,sys,time
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--map", default="${map}", help="")
parser.add_argument("--subset_map", default="${subset_map}", help="")
parser.add_argument("--out_map", default="${out_map}", help="")
args = parser.parse_args()

def subset_map(map, subset_map, out_map):
    '''
    Return: max and min base for each chromosome
    '''
    data = {}
    out = gzip.open(out_map, 'wb')

    for line in gzip.open(map):
        if 'chr' in line and 'position' in line:
            out.writelines(line)
        dat = line.strip().split()
        try:
            chr = dat[0]
            pos = dat[1]
            chr_pos = chr+"_"+pos
            data[chr_pos] = line
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
    subset_map(args.map, args.subset_map, args.out_map)

