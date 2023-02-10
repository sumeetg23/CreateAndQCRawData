from snakemake.utils import min_version, validate
import os
import sys
from helper import validateSamplesheet, validateOutput
import pandas as pd
import re
import collections
import json
import glob
import itertools

sys.path.append(os.path.abspath(os.getcwd()))
dirname=os.path.dirname
p = os.path.abspath(".")

min_version("6.4.0") # use-envmodules are not working in versions of 5.10.0 and below for clusters

NAMES ={}
SAMPLES =[]
TARGET = []
TARGET_BASE = []
LIBRARIES = []

#Snakemake configs and setup
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")
samplesheet = config["bcl2fastq"]["SampleSheet"]
Runfolder = config["bcl2fastq"]["OutputFolder"]
fastqre   = re.compile(r'\.fastq.gz$')

def validateBefore(outputfolder):
	success = validateSamplesheet(samplesheet)
	return True if success is not None else False

def getParentDir(wildcards):
	return dirname(str(config["bcl2fastq"]["SampleSheet"]))

def defineFileNames():
	lineno = 1
	for line in f:
		if 'Sample_ID' in line:
			column = line.split(",")
			sampleID    =column.index('Sample_ID')
			sampleName  =column.index('Sample_Name')
        else:
			column = line.split(",")
			SAMPLES += [column[sampleID]]
			NAMES[column[sampleID]] = column[sampleName]
			TARGET += [column[sampleID]+"_S"+str(lineno)+"_R1_001"]
			TARGET += [column[sampleID]+"_S"+str(lineno)+"_R2_001"]
			TARGET_BASE += [column[sampleID]+"_S"+str(lineno)]
	        lineno = lineno + 1

	IDS = tuple(TARGET)
	R1IDS = tuple(TARGET_BASE)
	
	return IDS, R1IDS

localrules: all

rule all:
    input:
        "multiqc_report.html"

rule bcl2fastq:
    input:
        input = config["bcl2fastq"]["SampleSheet"]
    params:
        barcode_mismatches = config["bcl2fastq"]["barcode_mismatch"],
        threads = config["bcl2fastq"]["threads"],
        infolder = config["bcl2fastq"]["RunFolder"],
        additionalOptions=[" "+config["bcl2fastq"]["other_params"],""][len(config["bcl2fastq"]["other_params"])>0],
    log:
        config["bcl2fastq"]["RunFolder"]+"/logs/e_bcl.log"
    output:
        bcl2fastqOutput = config["bcl2fastq"]["RunFolder"]+"/Stats/Stats.json"
		all = expand(config["bcl2fastq"]["RunFolder"]+"FASTQ/{sample}.fastq.gz", sample = IDS),
    threads: config["bcl2fastq"]["threads"]
    message:
        "Running bcl2fastq"
    shell:
        """
        bcl2fastq -R {params.infolder} --sample-sheet {input[0]} {params.additionalOptions} --barcode-mismatches {params.barcode_mismatches} -r {params.threads} -p {params.threads} >> {log} 2>&1
        """

rule fastqc:
    input:
        fq1 = "demultiplexed_reads/{sample}_R1_001.fastq.gz",
        fq2 = "demultiplexed_reads/{sample}_R2_001.fastq.gz"
    output:
        fq1 = "demultiplexed_reads/{sample}_R1_001_fastqc.html",
        fq2 = "demultiplexed_reads/{sample}_R2_001_fastqc.html"
    run:
        shell("""
        fastqc {input.fq1} 
        fastqc {input.fq2} 
        """)

rule multiqc:
    input:
        expand("demultiplexed_reads/{sample}_R1_001_fastqc.html", sample = R1IDS),
    output:
        report("multiqc_report.html")
    shell:
        "multiqc -f . ;"

onsuccess:
    print("Workflow finished without errors")

onerror:
    print("An error occurred during runtime")

onstart:
    print("Setting up and running pipeline")
