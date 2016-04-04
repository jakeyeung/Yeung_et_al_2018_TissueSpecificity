#!/usr/bin/env python
'''
DESCRIPTION

    Merge DHS bedfiles for easier R processing
    Example indir:
     /home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_stringent


FOR HELP

    python merge_dhs.py --help

AUTHOR:      Jake Yeung (jake.yeung@epfl.ch)
LAB:         Computational Systems Biology Lab (http://naef-lab.epfl.ch)
CREATED ON:  2016-04-04
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys
import argparse
import os
import csv


class BedRow():
    def __init__(self, row):
        '''
        '''
        self.chromo, self.start, self.end = row[0], row[1], row[2]
        self.val = row[3]
        self.coord = ':'.join([self.chromo, '-'.join([self.start, self.end])])


def coord_to_list(coord):
    '''
    chr1:123-345 -> [chr1, 123, 345]
    '''
    coord = coord.replace(":", "-")
    return(coord.split("-"))


def update_peak_dic(bed, tissue, dic, tissues):
    with open(bed, 'rb') as bedfile:
        bedreader = csv.reader(bedfile, delimiter='\t')
        for row in bedreader:
            Row = BedRow(row)
            if Row.coord not in dic:
                dic[Row.coord] = {}  # init for each tissue
                for tiss in tissues:
                    dic[Row.coord][tiss] = 0
            dic[Row.coord][tissue] = Row.val
    return(dic)


def write_dic_to_output(dic, coords, outf):
    with open(outf, 'wb') as outfile:
        outwriter = csv.writer(outfile, delimiter='\t')
        for coord, coord_dic in dic.iteritems():
            for tiss, val in coord_dic.iteritems():
                outrow = coord_to_list(coord)
                for vt in [val, tiss]:
                    outrow.append(vt)
                outwriter.writerow(outrow)


def main():
    parser = argparse.ArgumentParser(description='Merge DHS bedfiles for '
                                     'easier R processing')
    parser.add_argument('indir', metavar='INDIR',
                        help='indir')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='outfile')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    args = parser.parse_args()

    # store command line arguments for reproducibility
    CMD_INPUTS = ' '.join(['python'] + sys.argv)    # easy printing later
    # store argparse inputs for reproducibility / debugging purposes
    args_dic = vars(args)
    ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.iteritems()]
    ARG_INPUTS = ' '.join(ARG_INPUTS)

    # Print arguments supplied by user
    if not args.quiet:
        print('Command line inputs:')
        print(CMD_INPUTS)
        print ('Argparse variables:')
        print(ARG_INPUTS)

    beds = os.listdir(args.indir)
    beds = [b for b in beds if b.endswith(".bed")]
    tissues = [b.split(".")[0] for b in beds]
    paths = ['/'.join([args.indir, b]) for b in beds]
    print(beds)
    print(tissues)
    print(paths)

    peak_dic = {}
    for bedpath, tissue in zip(paths, tissues):
        dic = update_peak_dic(bedpath, tissue, peak_dic, tissues)

    coords = sorted(dic.keys())
    print(coords[0:5])
    # write dic to output
    write_dic_to_output(dic, coords, args.outfile)


if __name__ == '__main__':
    main()
