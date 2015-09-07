#!/usr/bin/env python
'''
DESCRIPTION

    After running dhs_merged_find_peaks.R (via run_dhs_merged_find_peaks.sh),
    merge the DHS counts together across all tissues into a single matrix.

FOR HELP

    python merge_dhs_counts.py --help

AUTHOR:      Jake Yeung (jake.yeung@epfl.ch)
LAB:         Computational Systems Biology Lab (http://naef-lab.epfl.ch)
CREATED ON:  2015-09-07
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse
import csv

def get_tissue_beds(bed_dirs_file):
    '''
    Input:
        Textfile with columns:
            1: Sample name (tissue)
            2: Path to bed file
    Output:
        Dic: {Tissue: Bedfile}
    '''
    beds_list = []
    with open(bed_dirs_file, 'rb') as readfile:
        jreader = csv.reader(readfile, delimiter = '\t')
        for row in jreader:
            print row
            colname, jpath = row[0], row[1]
            beds_list.append((colname, jpath))
    return beds_list

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('infile', metavar='INFILE',
                        help='File with two columns. First column contains\
                        column names (e.g., Liver) and second column contains\
                        path to filtered normalized scaled bed file')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='Matrix of DHS counts. First 3 columns are BED-style')
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

    tissue_beds_list = get_tissue_beds(args.infile)

    # open files
    infiles_list = []
    for (tiss, path) in tissue_beds_list:
        rfile = open(path, 'rb')
        reader_obj = csv.reader(rfile, delimiter = '\t')
        infiles_list.append((tiss, reader_obj))

    done = False
    # write to outfile
    with open(args.outfile, 'wb') as outfile:
        jwriter = csv.writer(outfile, delimiter = '\t')
        # write header
        header = ['chromo', 'start', 'end']  # init
        for (tiss, _) in tissue_beds_list:
            header.append(tiss)
        jwriter.writerow(header)
        # read files line by line
        writecount = 0
        while True:
            bedinforow = []
            countsrow = []
            for (_, readobj) in infiles_list:
                try:
                    row = readobj.next()
                except StopIteration:
                    done = True
                    break
                if len(bedinforow) == 0:
                    # add chromo start end to writerow only ONCE
                    bedinfo = row[0:3]
                    for bi in bedinfo:
                        bedinforow.append(bi)
                counts = float(row[3])  # bed format, 4th col is counts
                countsrow.append(counts)
            # check we are done, break out if so
            if done:
                break
            # only write if sum is not 0
            if sum(countsrow) > 0:
                jwriter.writerow(bedinforow + countsrow)
                writecount += 1
            else:
                pass

    print('%s rows written to: %s' % (writecount, args.outfile))

if __name__ == '__main__':
    main()
