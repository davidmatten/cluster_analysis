#!/usr/bin/python

import sys, os
import argparse


def main(infile, biomeasures, outfile):
# read files into 
    print infile
    print biomeasures
    print outfile
    print "end."



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Performs clustering analysis on tab delimited aligned protein sequences matching some biological measurement. This was developed using cd4 and log viral load. Example usage: python cluster_analysis.py -seqs='/path/to/sequence/file.txt' -bio_measure='/path/to/biological/measures.txt' -outf='/path/to/out.out'""")
    parser.add_argument('-seqs', '--sequences', type=str,
            help='The path to the tab delimited aligned protein sequence files. eg: "/path/to/sequences/files.txt"', required=True)
    parser.add_argument('-bio', '--biomeasures', type=str,
            help="The path to the tab delimited file containing biological measurements of the associated protein sequences in the sequences file.", required=True)
    parser.add_argument('-outf', '--outfile', type=str,
                        help = 'Path (including file name) to where outfile will be generated.', required=True)
    args = parser.parse_args()
    
    seqs = args.sequences
    biomeasures = args.biomeasures
    outfile = args.outfile

    if (seqs is None) or (biomeasures is None) or (outfile is None):
        print "Please specify all required arguments. infile, bio measures, outfile"
        sys.exit()

    if not os.path.exists(seqs):
        print "Please supply a valid path to the sequence source files"
        sys.exit()
    if not os.path.exists(biomeasures):
        print "Please supply a valid path to the bio measures source files"
        sys.exit()

    main(seqs, biomeasures, outfile)

