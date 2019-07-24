#!/usr/bin/env python
"""
Reads map file, chunk size
Returns file with chromosome chunk_start chunk_end
"""


def get_coordinates(bedFile, bed_output):
    '''
    Return: max and min base for each chromosome
    '''
    data = {}
    output = open(bed_output, 'w')
    for line in open(bedFile):
        line = line.strip().split()
        chrm = line[0]
        start = int(line[1])
        end = int(line[2])
        if chrm not in data:
            data[chrm] = [start, end]
        else:
            if start <= data[chrm][0]:
                data[chrm][0] = start
            if end > data[chrm][1]:
                data[chrm][1] = end
    for chrm in data:
        start = str(data[chrm][0])
        start_ = str(data[chrm][0])
        end = str(data[chrm][1])
        output.writelines(chrm+"\\t"+start_+"\\t"+end+"\\t"+chrm+" dna:chromosome chromosome:GRCh37:"+chrm+":"+start+":"+end+":1\\n")
    output.close()


if __name__ == '__main__':
    bedFile = "${bedFile}"
    bed_output = "${bed_output}"
    get_coordinates(bedFile, bed_output)
