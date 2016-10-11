'''
count_mismatches.py - count the number of mismatches per gene
====================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Count the number of high quality mismatches per gene per base sequenced. Will discard reads marked as duplicate

Usage
-----

.. Example use case

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
from collections import defaultdict

from CGAT import Experiment as E
from CGAT import GTF
from CGAT import IOTools

import pysam

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-b", "--bamfile", dest="bam", type="string",
                      help="BAM formated alignment file to test. Should have MD and NH tags set")
    parser.add_option("-t", "--quality-threshold", dest="threshold", type="int",
                       default=30,
                       help="minimum quality threshold for a mismatched base to count")


    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    bamfile = pysam.AlignmentFile(options.bam)
    options.stdout.write("\t".join(["gene_id",
                                    "mismatches",
                                    "bases",
                                    "low_qual",
                                    "a","t","c","g"]) + "\n")

    for gene in GTF.flat_gene_iterator(GTF.iterator(options.stdin)):

        start = min(e.start for e in gene)
        end = max(e.end for e in gene)

        reads = bamfile.fetch(gene[0].contig, start, end)

        gene_id = gene[0].gene_id
        mm_count = 0
        base_count = 0

        skipped = 0
        mm_bases = defaultdict(int)
        for read in reads:

            if read.is_unmapped:
                continue
            if read.is_duplicate:
                continue

            if read.get_tag("NH") > 1:
                continue

            alignment = read.get_aligned_pairs(matches_only=True, with_seq=True)
            base_count += sum(1 for b in read.query_alignment_qualities
                              if b >= options.threshold)
            if read.get_tag("NM") == 0:
                continue

            # mismatches
            mismatches = [base for base in alignment
                          if base[2].islower() and
                          start <= base[1] < end]

            qualities = read.query_qualities
            
            hq_mm = sum(1 for base in mismatches
                        if qualities[base[0]] >= options.threshold)

            for base in mismatches:
                mm_bases[base[2]] += 1

            mm_count += hq_mm
            skipped += len(mismatches) - hq_mm

        outline = "\t".join(map(str,[gene_id,
                                     mm_count,
                                     base_count,
                                     skipped,
                                     mm_bases['a'],
                                     mm_bases['t'],
                                     mm_bases['c'],
                                     mm_bases['g']]))
        options.stdout.write(outline + "\n")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
