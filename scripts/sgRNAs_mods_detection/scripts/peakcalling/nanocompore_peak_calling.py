#!/usr/bin/env python

import argparse, sys, collections
import numpy as np
import re
from scipy.signal import find_peaks

###############################################################################
def parse_args(args):
    '''
    parses command line arguments
    '''
    parser = argparse.ArgumentParser(description ='Filters bedtools intersect output file')

    parser.add_argument('--infile', '-i', type = str, help = 'input bed intersection file')

    parser.add_argument('--outfile', '-o', type = str, help = 'output bed file')

    return parser.parse_args()
###############################################################################

###############################################################################
def convert_nanocompore_tsv_to_bed(options):
    '''
    '''
    infile = options.infile
    outfile = options.outfile
    with open(outfile, 'w') as out:
        for orf, pos, gpos in peakcall_nanocompore_transcript_pvalues(infile):
            out.write(f"{orf}\t{gpos}\t{gpos}\t.\t.\t+\t{pos}\n")
###############################################################################

###############################################################################
def peakcall_nanocompore_transcript_pvalues(infile):
    '''
    '''
    for orf, gmm_pvalues, positions, gpositions, kmers, lors in collect_transcript_information(infile):
        gmm_pvalues, positions, gpositions, kmers, lors = filter_for_u(gmm_pvalues, positions, gpositions, kmers, lors)
        gmm_pvalues = -1*np.log10(gmm_pvalues)
        passing_gmm_pvalues = []
        peak_idx = find_peaks(gmm_pvalues, height=2, distance=5)

        for idx in peak_idx[0]:
            if abs(lors[idx]) >= 0.5 and gmm_pvalues[idx] >= 2:
                yield orf, positions[idx], gpositions[idx]
###############################################################################

###############################################################################
def filter_for_u(gmm_pvalues, positions, gpositions, kmers, lors):
    '''
    '''
    gmm_pvalues_of_Us = []
    positions_of_Us = []
    gpositions_of_Us = []
    kmers_of_Us = []
    lors_of_Us = []
    for gmm, pos, gpos, kmer, lor in zip(gmm_pvalues, positions, gpositions, kmers, lors):
        if 'U' in kmer or 'T' in kmer:
            gmm_pvalues_of_Us.append(gmm)
            positions_of_Us.append(pos)
            gpositions_of_Us.append(gpos)
            kmers_of_Us.append(kmer)
            lors_of_Us.append(lor)

    return gmm_pvalues_of_Us, positions_of_Us, gpositions_of_Us, kmers_of_Us, lors_of_Us
###############################################################################

###############################################################################
def collect_transcript_information(infile):
    '''
    '''
    current_orf = ''
    gmm_pvalues = []
    positions = []
    gpositions = []
    kmers = []
    lors = []
    for orf, gmm_pvalue, lor, pos, gpos, kmer in read_nanocompore_tsv(infile):
        if current_orf and orf != current_orf:
            yield current_orf, gmm_pvalues, positions, gpositions, kmers, lors
            current_orf = orf
            gmm_pvalues = [gmm_pvalue]
            positions = [pos]
            gpositions = [gpos]
            kmers = [kmer]
            lors = [lor]
        else:
            current_orf = orf
            gmm_pvalues.append(gmm_pvalue)
            positions.append(pos)
            gpositions.append(gpos)
            kmers.append(kmer)
            lors.append(lor)

    yield current_orf, gmm_pvalues, positions, gpositions, kmers, lors
###############################################################################

###############################################################################
def read_nanocompore_tsv(infile):
    '''
    '''
    with open(infile, 'r') as tsv:
        line = tsv.readline()
        for line in tsv:
            line = line.strip().split('\t')
            try:
                gmm_pvalue = float(line[6])
            except:
                gmm_pvalue = float('NaN')

            try:
                lor = float(line[-1])
            except:
                lor = float('NaN')

            orf = line[3].strip()
            pos = line[0].strip()
            gpos = line[2].strip()
            kmer = line[5].strip()
            yield orf, gmm_pvalue, lor, pos, gpos, kmer
###############################################################################

###############################################################################
def main(args):
    #Parse the inputs args/options
    options = parse_args(args)

    convert_nanocompore_tsv_to_bed(options)
###############################################################################

if (__name__ == "__main__"):
    main(sys.argv)
    raise SystemExit

