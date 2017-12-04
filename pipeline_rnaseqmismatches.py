##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################
"""===========================
Pipeline template
===========================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

.. Replace the documentation below with your own description of the
   pipeline's purpose

Overview
========

This pipeline computes the word frequencies in the configuration
files :file:``pipeline.ini` and :file:`conf.py`.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_rnaseqmismatches.py config

Input files
-----------

None required except the pipeline configuration files.

Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

.. Add any additional external requirements such as 3rd party software
   or R modules below:

Requirements:

* samtools >= 1.1

Pipeline output
===============

.. Describe output files of the pipeline here

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys
import os
import sqlite3
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import re

# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS["projectsrc"] = os.path.dirname(__file__)
#for key, value in PARAMS.iteritems():
#    print "%s:\t%s" % (key,value)


# add configuration values from associated pipelines
#
# 1. pipeline_annotations: any parameters will be added with the
#    prefix "annotations_". The interface will be updated with
#    "annotations_dir" to point to the absolute path names.
PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__",
    prefix="annotations_",
    update_interface=True))


# if necessary, update the PARAMS dictionary in any modules file.
# e.g.:
#
# import CGATPipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS
#
# Note that this is a hack and deprecated, better pass all
# parameters that are needed by a function explicitely.

# -----------------------------------------------
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


# ---------------------------------------------------
# Specific pipeline tasks

@follows(mkdir("readgroups.dir"))
@transform("*.bam",regex(r"(.+).bam"),r"readgroups.dir/\1.readgroups.bam")
def add_read_groups(infile, outfile):
    platform = PARAMS["platform"]
    groupsample = PARAMS["groupsample"]
    statement = '''java -jar /shared/sudlab1/General/apps/bio/picard-tools-1.135/picard.jar
                   AddOrReplaceReadGroups
                   I=%(infile)s
                   O=%(outfile)s
                   RGLB=lib1
                   RGPL=%(platform)s
                   RGPU=unit1
                   RGSM=%(groupsample)s'''

    job_memory = "4G"
    P.run()



@follows(mkdir("deduped.dir"))
@transform(add_read_groups,
           regex(r"readgroups.dir/(.+).readgroups.bam"),
           r"deduped.dir/\1.bam")
def dedup_bams(infile, outfile):
    '''Use MarkDuplicates to mark dupliceate reads'''
    tempfile=P.snip(outfile, ".bam") + ".temp.bam"   
    metrics=P.snip(outfile, ".bam") + ".metrics.tsv"
    temporary = PARAMS["tmpdir"]
    statement = '''MarkDuplicates I=%(infile)s
                                  O=%(tempfile)s
                                  M=%(metrics)s
                                  TMP_DIR=%(temporary)s > %(outfile)s.log;

                    checkpoint;
                                samtools view 
                                -F 1024
                                -b
                                %(tempfile)s
                                > %(outfile)s;
                  
                    checkpoint;
                                rm -r %(tempfile)s

                    checkpoint;

                                samtools index %(outfile)s'''

    job_memory = "15G"
    P.run()

@active_if(not(PARAMS['vcfavail']))
@follows(mkdir("split.dir"))
@transform(dedup_bams,regex(r"deduped.dir/(.+).bam"),r"split.dir/\1.split.bam")
def splitbams(infile,outfile):
    '''use GATK splitNcigar to split reads into exon segements'''
    fasta = os.path.join(PARAMS["fasta"],PARAMS["genome"]) + ".fasta"
    fastamap = PARAMS["mapfasta"]
    drctry= PARAMS["tmpdir"]
    statement = '''java -Xmx10G -Djava.io.tmpdir=%(drctry)s -jar ~/Downloads/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
                   -T SplitNCigarReads 
                   -R %(fastamap)s
                   -I %(infile)s 
                   -o %(outfile)s 
                   -rf ReassignOneMappingQuality 
                   -RMQF 255 
                   -RMQT 60 
                   -U ALLOW_N_CIGAR_READS'''

    job_memory = "12G"
    P.run()

#@follows(mkdir("BaseRecalibration.dir"))
#@transform(splitbams,regex(r"split.dir/(.+).split.bam"),r"BaseRecalibration.dir/\1.recal.bam")
#def baserecal(infile,outfile):
#    fasta = os.path.join(PARAMS["fasta"],PARAMS["genome"]) + ".fasta"
#    fastamap = PARAMS["mapfasta"]
#    drctry= PARAMS["tmpdir"]
#    statement = '''java -Xmx10G -Djava.io.tmpdir=%(drctry)s -jar ~/Downloads/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar 
#                   -T BaseRecalibrator
#                   -R %(fastamap)s
#                   -I %(infile)s 
#                   -o %(outfile)s
#                   '''
#  
#    job_memory = "12G"
#    P.run()               

@follows(mkdir("Variantcalls.dir"))
@transform(splitbams,regex(r"split.dir/(.+).split.bam"),r"Variantcalls.dir/\1.vcf.gz")
def variantcalling(infile,outfile):
    fasta = os.path.join(PARAMS["fasta"],PARAMS["genome"]) + ".fasta"
    fastamap = PARAMS["mapfasta"]
    drctry= PARAMS["tmpdir"]
    tempfile=P.snip(outfile,".gz")
    statement = '''java -Xmx10G -Djava.io.tmpdir=%(drctry)s -jar ~/Downloads/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar 
                   -T HaplotypeCaller
                   -R %(fastamap)s 
                   -I %(infile)s 
                   -dontUseSoftClippedBases 
                   -stand_call_conf 20.0 
                   -o %(tempfile)s;
                   checkpoint;
                   bgzip %(tempfile)s;
                   '''  

    job_memory = "12G"
    P.run()                 


@transform(variantcalling,regex(r"(.+).vcf.gz"),r"\1.reheader.vcf.gz.tbi")
def renamesample(infile,outfile):
    renamefile=P.snip(infile,".vcf.gz")+".txt"
    tempfile=P.snip(outfile,".tbi")
    vcfpattern=PARAMS["vcfpattern"]
    statement='''touch %(renamefile)s;
                 checkpoint;
                 echo $(echo %(infile)s | egrep -o %(vcfpattern)s) > %(renamefile)s;
                 checkpoint;
                 bcftools reheader -s %(renamefile)s %(infile)s > %(tempfile)s;
                 checkpoint;
                 tabix -p vcf %(tempfile)s'''

    job_memory = "4G"
    P.run()

#add_inputs("geneset_all.gtf.gz")
# ---------------------------------------------------
@active_if(PARAMS['vcfavail'])
@follows(mkdir("mismatches.dir"))
@transform(dedup_bams,
           regex(r"deduped.dir/(.+).bam"),
           r"mismatches.dir/\1.tsv.gz")
def count_mismatches(infile, outfile):
    ''' Count mismatches per sequenced base, per read, discarding duplicated reads
    and low quality bases'''
    fastapath = os.path.join(PARAMS["fasta"],PARAMS["genome"])
    vcfpath = PARAMS["vcf"]
    gtfpath = PARAMS["gtf"]
    redipath = PARAMS["redipath"]
    sampat = PARAMS["samplepattern"]
    samplepattern = '"%s"'%(sampat)
    quality_threshold = PARAMS["quality_threshold"]
    statement = '''python %(projectsrc)s/count_mismatches.py
                                         -I %(gtfpath)s
                                         --bamfile=%(infile)s
                                         --quality-threshold=%(quality_threshold)s
                                         --fasta-path=%(fastapath)s
                                         -d %(samplepattern)s
                                         --vcf-path=%(vcfpath)s
                                         --REDI-path=%(redipath)s
                                         -S %(outfile)s
                                         -L %(outfile)s.log
                                         -v5 '''
    job_memory="6G"
    P.run()
#@transform(dedup_bams,
#           formatter(),
#           "mismatches.dir/{basename[0]}.tsv.gz")

@active_if(not(PARAMS['vcfavail']))
@follows("renamesample")
@follows(mkdir("mismatches.dir"))
@transform(dedup_bams,
           regex(r"deduped.dir/(.+).bam"),
           r"mismatches.dir/\1.tsv.gz")
def count_mismatches_with_VCF(infile, outfile):
    ''' Count mismatches per sequenced base, per read, discarding duplicated reads
    and low quality bases'''
    fastapath = os.path.join(PARAMS["fasta"],PARAMS["genome"])
    vcfname = re.search(r"deduped.dir/(.+).bam", infile, flags = 0).group(1) + ".reheader.vcf.gz"
    vcfpath = "Variantcalls.dir/" + vcfname
    gtfpath = PARAMS["gtf"]
    redipath = PARAMS["redipath"]
    sampat = PARAMS["samplepattern"]
    samplepattern = '"%s"'%(sampat)
    quality_threshold = PARAMS["quality_threshold"]
    statement = '''python %(projectsrc)s/count_mismatches.py
                                         -I %(gtfpath)s
                                         --bamfile=%(infile)s
                                         --quality-threshold=%(quality_threshold)s
                                         --fasta-path=%(fastapath)s
                                         --vcf-path=%(vcfpath)s
                                         --REDI-path=%(redipath)s
                                         -d %(samplepattern)s
                                         -S %(outfile)s
                                         -L %(outfile)s.log
                                         -v5'''
    job_memory="8G"
    P.run()
# ---------------------------------------------------
@merge([count_mismatches,count_mismatches_with_VCF],
       "mismatch_counts.load")
def merge_mismatch_counts(infiles, outfile):
    '''Load the results of mismatch counting into the database'''
    
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="mismatches.dir/(.+)-(.+).tsv.gz",
                         cat="tissue,replicate",
                         options="-i tissue -i replicate -i gene_id",
			 job_memory="10G")


# ---------------------------------------------------
# Generic pipeline tasks
@follows(merge_mismatch_counts)
def full():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.

    Any existing report will be overwritten.
    '''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.

    This will update a report with any changes inside the report
    document or code. Note that updates to the data will not cause
    relevant sections to be updated. Use the cgatreport-clean utility
    first.
    '''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report in the CGAT downloads directory.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
