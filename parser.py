#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Parses eggNOG text files and the supplied BLAST result XML file to produce
summary statistics about the matches to eggNOG clusters.

Sample command line to obtain an XML file of hits to the bactNOG collection
of sequences:
blastall -a 16 -b 1 -v 1 -p blastp -e 0.001 -m 7 -d sequences/all.fa -i genome.faa -o genome.results.xml
-a 16: number of threads
-b 1: number of alignments to show
-v 1: number of descriptions to show
-p: program to use
-e: expectation threshold
-m 7: use XML format for output
-d PATH: to the database containing all bactNOG sequences (merged manually; must be format-ted first)
-i genome: analyzed genome, multifasta of proteins
-o filename: where to write the results to
'''


import sys
from argparse import ArgumentParser


# Mapping of single-letter category names to their full names.
oneletter_to_full = {'J' : 'Translation, ribosomal structure and biogenesis',
                     'A' : 'RNA processing and modification',
                     'K' : 'Transcription',
                     'L' : 'Replication, recombination and repair',
                     'B' : 'Chromatin structure and dynamics',
                     'D' : 'Cell cycle control, cell division, chromosome partitioning',
                     'Y' : 'Nuclear structure',
                     'V' : 'Defense mechanisms',
                     'T' : 'Signal transduction mechanisms',
                     'M' : 'Cell wall/membrane/envelope biogenesis',
                     'N' : 'Cell motility',
                     'Z' : 'Cytoskeleton',
                     'W' : 'Extracellular structures',
                     'U' : 'Intracellular trafficking, secretion, and vesicular transport',
                     'O' : 'Posttranslational modification, protein turnover, chaperones',
                     'C' : 'Energy production and conversion',
                     'G' : 'Carbohydrate transport and metabolism',
                     'E' : 'Amino acid transport and metabolism',
                     'F' : 'Nucleotide transport and metabolism',
                     'H' : 'Coenzyme transport and metabolism',
                     'I' : 'Lipid transport and metabolism',
                     'P' : 'Inorganic ion transport and metabolism',
                     'Q' : 'Secondary metabolites biosynthesis, transport and catabolism',
                     'R' : 'General function prediction only',
                     'S' : 'Function unknown'
                     }

# Mapping of single-letter category names to super-category names.
oneletter_to_super = {'J': 'INFORMATION STORAGE AND PROCESSING',
                      'A': 'INFORMATION STORAGE AND PROCESSING',
                      'K': 'INFORMATION STORAGE AND PROCESSING',
                      'L': 'INFORMATION STORAGE AND PROCESSING',
                      'B': 'INFORMATION STORAGE AND PROCESSING',
                      'D': 'CELLULAR PROCESSES AND SIGNALING',
                      'Y': 'CELLULAR PROCESSES AND SIGNALING',
                      'V': 'CELLULAR PROCESSES AND SIGNALING',
                      'T': 'CELLULAR PROCESSES AND SIGNALING',
                      'M': 'CELLULAR PROCESSES AND SIGNALING',
                      'N': 'CELLULAR PROCESSES AND SIGNALING',
                      'Z': 'CELLULAR PROCESSES AND SIGNALING',
                      'W': 'CELLULAR PROCESSES AND SIGNALING',
                      'U': 'CELLULAR PROCESSES AND SIGNALING',
                      'O': 'CELLULAR PROCESSES AND SIGNALING',
                      'C': 'METABOLISM', 'G': 'METABOLISM', 'E': 'METABOLISM',
                      'F': 'METABOLISM', 'H': 'METABOLISM', 'I': 'METABOLISM',
                      'P': 'METABOLISM', 'Q': 'METABOLISM',
                      'R': 'POORLY CHARACTERIZED', 'S': 'POORLY CHARACTERIZED'
                      }

# Format string for the final statistics display.
final = '''
INFORMATION STORAGE AND PROCESSING: %s
 [J] Translation, ribosomal structure and biogenesis: %s 
 [A] RNA processing and modification: %s
 [K] Transcription: %s
 [L] Replication, recombination and repair: %s
 [B] Chromatin structure and dynamics: %s

CELLULAR PROCESSES AND SIGNALING: %s
 [D] Cell cycle control, cell division, chromosome partitioning: %s
 [Y] Nuclear structure: %s
 [V] Defense mechanisms: %s
 [T] Signal transduction mechanisms: %s
 [M] Cell wall/membrane/envelope biogenesis: %s
 [N] Cell motility: %s
 [Z] Cytoskeleton: %s
 [W] Extracellular structures: %s
 [U] Intracellular trafficking, secretion, and vesicular transport: %s
 [O] Posttranslational modification, protein turnover, chaperones: %s

METABOLISM: %s
 [C] Energy production and conversion: %s
 [G] Carbohydrate transport and metabolism: %s
 [E] Amino acid transport and metabolism: %s
 [F] Nucleotide transport and metabolism: %s
 [H] Coenzyme transport and metabolism: %s
 [I] Lipid transport and metabolism: %s
 [P] Inorganic ion transport and metabolism: %s
 [Q] Secondary metabolites biosynthesis, transport and catabolism: %s

POORLY CHARACTERIZED: %s
 [R] General function prediction only: %s
 [S] Function unknown: %s
'''

# (Empty) dict of per-category counters.
category_counters = {'J':0, 'A':0, 'K':0, 'L':0, 'B':0, 'D':0, 'Y':0, 'V':0, 'T':0, 'M':0,
                     'N':0, 'Z':0, 'W':0, 'U':0, 'O':0, 'C':0, 'G':0, 'E':0, 'F':0, 'H':0,
                     'I':0, 'P':0, 'Q':0, 'R':0, 'S':0 }


def parse_bactnog_members(fname):
    '''
    Given 'fname', the path to bactNOG.members.txt file, parses it into a dictionary
    where bactNOG sequence names are keys and bactNOG group numbers are values.
    '''
    # Sample lines: group number, sequence name, unknown, unknown
    #bactNOG00016    31964.CMS_2765    3    372
    #bactNOG00016    36870.WGLp486    1    361
    #bactNOG00016    41514.SARI_01986    1    362
    mapping = {}
    with open(fname) as fh:
        for line in fh:
            number, sequence = line.split(maxsplit = 2)
            mapping[sequence] = number
    return mapping


def parse_bactnog_funccat(fname):
    '''
    Given 'fname', the path to the bactNOG.funccat.txt file, return a dictionary
    where bactNOG group numbers are keys and single-letter category names are values.
    '''
    # Sample lines: group number, 1-letter category name
    #bactNOG59442    L
    #bactNOG41622    T
    #bactNOG41621    S
    mapping = {}
    with open(fname) as fh:
        for line in fh:
            number, category = line.strip().split()
            mapping[number] = category
    return mapping


def main(argv=None):
    '''Command line options.'''
    
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    # Setup argument parser
    parser = ArgumentParser()
    parser.add_argument(dest="paths", help="paths to folder(s) with source file(s) [default: %(default)s]", metavar="path", nargs='+')
    
    # Process arguments
    args = parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())
