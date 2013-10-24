#!/usr/bin/python
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
from Bio import SearchIO as SIO


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

information_codes_list = ['J', 'A', 'K', 'L', 'B']
signaling_codes_list = ['D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O']
metabolism_codes_list = ['C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q']
unknown_codes_list = ['R', 'S']

# Format string for the final statistics display.
final = '''
INFORMATION STORAGE AND PROCESSING: {s[INFORMATION STORAGE AND PROCESSING]}
 [J] Translation, ribosomal structure and biogenesis: {c[J]}
 [A] RNA processing and modification: {c[A]}
 [K] Transcription: {c[K]}
 [L] Replication, recombination and repair: {c[L]}
 [B] Chromatin structure and dynamics: {c[B]}

CELLULAR PROCESSES AND SIGNALING: {s[CELLULAR PROCESSES AND SIGNALING]}
 [D] Cell cycle control, cell division, chromosome partitioning: {c[D]}
 [Y] Nuclear structure: {c[Y]}
 [V] Defense mechanisms: {c[V]}
 [T] Signal transduction mechanisms: {c[T]}
 [M] Cell wall/membrane/envelope biogenesis: {c[M]}
 [N] Cell motility: {c[N]}
 [Z] Cytoskeleton: {c[Z]}
 [W] Extracellular structures: {c[W]}
 [U] Intracellular trafficking, secretion, and vesicular transport: {c[U]}
 [O] Posttranslational modification, protein turnover, chaperones: {c[O]}

METABOLISM: {s[METABOLISM]}
 [C] Energy production and conversion: {c[C]}
 [G] Carbohydrate transport and metabolism: {c[G]}
 [E] Amino acid transport and metabolism: {c[E]}
 [F] Nucleotide transport and metabolism: {c[F]}
 [H] Coenzyme transport and metabolism: {c[H]}
 [I] Lipid transport and metabolism: {c[I]}
 [P] Inorganic ion transport and metabolism: {c[P]}
 [Q] Secondary metabolites biosynthesis, transport and catabolism: {c[Q]}

POORLY CHARACTERIZED: {s[POORLY CHARACTERIZED]}
 [R] General function prediction only: {c[R]}
 [S] Function unknown: {c[S]}
'''

# (Empty) dict of per-category counters.
category_counters = {'J':0, 'A':0, 'K':0, 'L':0, 'B':0, 'D':0, 'Y':0, 'V':0, 'T':0, 'M':0,
                     'N':0, 'Z':0, 'W':0, 'U':0, 'O':0, 'C':0, 'G':0, 'E':0, 'F':0, 'H':0,
                     'I':0, 'P':0, 'Q':0, 'R':0, 'S':0 }

# (Empty) dict for super-category counters.
super_counters = {'CELLULAR PROCESSES AND SIGNALING': 0, 'POORLY CHARACTERIZED': 0,
                  'INFORMATION STORAGE AND PROCESSING': 0,  'METABOLISM': 0}


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
            number, sequence, unk1, unk2 = line.split()
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


def sequence2category(sequence_id, s2b, b2c):
    '''
    Given a sequence_id of the bactNOG database hit, returns one-letter bactNOG category.
    Returns None if some of the mappings are absent.
    s2b: sequence id to bactnog dictionary
    b2c: bactnog number to category dictionary
    '''
    if sequence_id in s2b:
        if s2b[sequence_id] in b2c:
            return b2c[s2b[sequence_id]]
    return None


def main():
    '''Command line options.'''
    # Setup argument parser
    parser = ArgumentParser()
    parser.add_argument(dest="blast_xml", metavar="BLAST_result.xml",
                        help="path to the BLAST XML result file")
    parser.add_argument(dest="funccat", metavar="bactNOG.funccat.txt", default='bactNOG.funccat.txt',
                        help="path to the bactNOG.funccat.txt file [default: %(default)s]")
    parser.add_argument(dest="members", metavar="bactNOG.members.txt", default='bactNOG.members.txt',
                        help="path to the bactNOG.members.txt file [default: %(default)s]")
    # Process arguments
    args = parser.parse_args()

    bactnog2category = parse_bactnog_funccat(args.funccat)
    seqid2bactnog = parse_bactnog_members(args.members)

    # ORFs which did not have a match.
    nohits_counter = 0
    # ORFs which had a match but were not assigned to categories.
    nocat_counter = 0
    # Header for hits.
    print "Query\tHit\tCategory\tSuper-category"
    for blast_result in SIO.parse(args.blast_xml, 'blast-xml'):
        seqid = blast_result.id.split('.')[0][8:]
        if len(blast_result) == 0:
            nohits_counter += 1
            print "%s\tNo hits" % seqid
            continue
        category = sequence2category(blast_result[0].id, seqid2bactnog, bactnog2category)
        if category != None:
            for i in range(len(category)):
                cat = category[i]
                category_counters[cat] += 1
                # Print a hit with category and super-category.
                print "%s\t%s\t%s\t%s" % (seqid, blast_result[0].id,
                                                oneletter_to_full[cat],
                                                oneletter_to_super[cat])
        else:
            nocat_counter += 1
            print "%s\t%s\tNo category" % (seqid, blast_result[0].id)

    # Sum up individual categories to get super-level stats.
    for code in information_codes_list:
        super_counters['INFORMATION STORAGE AND PROCESSING'] += category_counters[code]
    for code in signaling_codes_list:
        super_counters['CELLULAR PROCESSES AND SIGNALING'] += category_counters[code]
    for code in metabolism_codes_list:
        super_counters['METABOLISM'] += category_counters[code]
    for code in unknown_codes_list:
        super_counters['POORLY CHARACTERIZED'] += category_counters[code]

    print 'Queries with no hits:', nohits_counter
    print 'Queries with no categories:', nocat_counter
    print final.format(c=category_counters, s=super_counters)


if __name__ == "__main__":
    sys.exit(main())
