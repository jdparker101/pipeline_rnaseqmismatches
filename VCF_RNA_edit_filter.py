'''
VCF_RNA_edit_filter.py - Filter out RNA editing events from VCFs
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
import CGAT.Bed as Bed
import textwrap
import pysam
import vcf
import re


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])


    parser.add_option("-p", "--vcf-path", dest="vcfpath", type="string",
                       help="path to indexed vcf file for dataset  of choice")
   

    (options, args) = E.Start(parser, argv=argv)

    vcffile = vcf.Reader(open(options.vcfpath,"r"))

    for record in vcffile:
        if (record.REF == "A" and record.ALT[0] == "G") or (record.REF == "T" and record.ALT[0] == "C"):
            continue
        else:
            line=(str(record.REF),str(record.ALT[0]),str(record.POS))
            outline="\t".join(line)            
            options.stdout.write(outline + "\n")

if __name__ == "__main__":
    sys.exit(main(sys.argv))
            
        

