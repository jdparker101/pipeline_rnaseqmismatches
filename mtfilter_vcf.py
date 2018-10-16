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
    parser.add_option("-o", "--out", dest="out",
                      help="name of output file")
    (options, args) = E.Start(parser, argv=argv)

    vcf_writer = vcf.Writer(open('/dev/null', 'w'), vcf_reader)
    vcffile = vcf.Reader(open(options.vcfpath,"r"))
    vcfregion = vcffile.fetch("chrM")
    vcf_writer = vcf.Writer(open('%s'%(out),'w'), vcfregion)
    for record in vcfregion:
    vcf_writer.write_record(record)
    vcf_writer.flush()


