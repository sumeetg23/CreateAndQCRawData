#Author:	Sumeet Gupta
#Email:		sumeetg23@yahoo.com
# This is a pipeline to convert bcl files to fastq.
# 

from snakemake.utils import min_version, validate
import os
import sys

sys.path.append(os.path.abspath(os.getcwd())+"/scripts/")

from helper import validateSamplesheet, validateOutput
import pandas as pd
import re
import collections
import json
import glob
import itertools

import warnings
warnings.filterwarnings('ignore')

dirname=os.path.dirname
p = os.path.abspath(".")


min_version("5.19.3") # use-envmodules are not working in versions of 5.10.0 and below for clusters

NAMES ={}
LIBRARIES = []

#Snakemake configs and setup
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")
samplesheet = config["bcl2fastq"]["SampleSheet"]
Runfolder = config["bcl2fastq"]["RunFolder"]
fastqre   = re.compile(r'\.fastq.gz$')

def validateBefore(outputfolder):
	success = validateSamplesheet(samplesheet)
	return True if success is not None else False

def getParentDir(wildcards):
	return dirname(str(config["bcl2fastq"]["SampleSheet"]))

def defineFileNames():
	f = open(config["bcl2fastq"]["SampleSheet"], 'r')
	SAMPLES =[]
	LANES = []
	TARGET = []
	TARGET_BASE = []
	lineno=[]
	lineno = [0 for i in range(8)] 
	next(f)
	for line in f:
		if 'Sample_ID' in line:
			column = line.split(",")
			sampleID = column.index('Sample_ID')
			sampleName = column.index('Sample_Name')
			lanenum = column.index('Lane')
		else:
			column = line.split(",")
			lineno[int(column[lanenum])] = lineno[int(column[lanenum])] + 1
			SAMPLES += [column[sampleID]]
			TARGET += [column[sampleID]+"_S"+str(lineno[int(column[lanenum])])+"_L00"+column[lanenum]+"_R1_001"]
			TARGET += [column[sampleID]+"_S"+str(lineno[int(column[lanenum])])+"_L00"+column[lanenum]+"_R2_001"]
			TARGET_BASE += [column[sampleID]+"_S"+str(lineno[int(column[lanenum])])+"_L00"+column[lanenum]]
			LANES.append(column[lanenum])

	return tuple(TARGET), tuple(TARGET_BASE), tuple(LANES)

IDS, R1IDS, LC = defineFileNames()

LC = set(LC)

print(R1IDS)

localrules: all, multiqc

rule all:
	input:
		config["bcl2fastq"]["RunFolder"]+"Stats/AdapterTrimming.txt",
		expand(config["bcl2fastq"]["RunFolder"]+"FASTQC/{sample}_R1_001_fastqc.html", sample = R1IDS),
		expand(config["bcl2fastq"]["RunFolder"]+"FASTQSCREEN/{sample}_R1_001_screen.txt", sample = R1IDS),
		expand(config["bcl2fastq"]["RunFolder"]+"FASTQSCREEN/{sample}_R1_001_screen.png", sample = R1IDS),
		expand(config["bcl2fastq"]["RunFolder"]+"MULTIQC/FASTQC_MULTIQC_L{lanenumber}_Interactive.html", lanenumber = LC),
		expand(config["bcl2fastq"]["RunFolder"]+"MULTIQC/FASTQSCREEN_MULTIQC_L{lanenumber}_Interactive.html", lanenumber = LC),
		Runfolder+"MULTIQC/multiqc_report.html"

rule bcl2fastq:
	input:
		input = config["bcl2fastq"]["SampleSheet"]
	params:
		barcode_mismatches = config["bcl2fastq"]["barcode_mismatch"],
		threads = config["bcl2fastq"]["threads"],
		infolder = config["bcl2fastq"]["RunFolder"],
		additionalOptions=[" "+config["bcl2fastq"]["other_params"],""][len(config["bcl2fastq"]["other_params"])>0]
	log:
		config["bcl2fastq"]["RunFolder"]+"logs/e_bcl.log"
	output:
		config["bcl2fastq"]["RunFolder"]+"Stats/AdapterTrimming.txt",
		bcl2fastqfiles = expand(config["bcl2fastq"]["RunFolder"]+"FASTQ/{sample}_R1_001.fastq.gz", sample = R1IDS)
	threads: config["bcl2fastq"]["threads"]
	message:
		"Running bcl2fastq"
	shell:
		"bcl2fastq -R {params.infolder} --output-dir={params.infolder} --sample-sheet {input[0]} {params.additionalOptions} --barcode-mismatches {params.barcode_mismatches} >> {log}"

rule fastqc:
	input:
		rt = Runfolder+"Stats/AdapterTrimming.txt",
		fq1 = config["bcl2fastq"]["RunFolder"]+"FASTQ/{sample}_R1_001.fastq.gz"
	output:
		config["bcl2fastq"]["RunFolder"]+"FASTQC/{sample}_R1_001_fastqc.html"
	shell:
		"fastqc --noextract --format fastq --threads 4 -o {Runfolder}FASTQC/ {input.fq1}"

rule fastqscreen:
	input:
		rt = Runfolder+"Stats/AdapterTrimming.txt",
		FASTQSCREENdir = Runfolder,
		fq1 = config["bcl2fastq"]["RunFolder"]+"FASTQ/{sample}_R1_001.fastq.gz"
	output:
		config["bcl2fastq"]["RunFolder"]+"FASTQSCREEN/{sample}_R1_001_screen.txt",
		config["bcl2fastq"]["RunFolder"]+"FASTQSCREEN/{sample}_R1_001_screen.png"
	shell:
		"fastq_screen --conf /usr/share/FastQ-Screen-0.15.2/fastq_screen.conf --outdir {input.FASTQSCREENdir}/FASTQSCREEN/ --aligner bowtie2 --force --subset 1000000 {input.fq1}"

rule multiqc:
	input:
		expand(config["bcl2fastq"]["RunFolder"]+"FASTQC/{sample}_R1_001_fastqc.html", sample = R1IDS),
		expand(config["bcl2fastq"]["RunFolder"]+"FASTQSCREEN/{sample}_R1_001_screen.png", sample = R1IDS),
		run = Runfolder
	output:
		expand(config["bcl2fastq"]["RunFolder"]+"MULTIQC/FASTQC_MULTIQC_L{lanenumber}_Interactive.html", lanenumber = LC),
		expand(config["bcl2fastq"]["RunFolder"]+"MULTIQC/FASTQSCREEN_MULTIQC_L{lanenumber}_Interactive.html", lanenumber = LC),
		report(Runfolder+"FASTQC/multiqc_report.html")
	shell:
		"PreProcess-Summary.sh -r {input.run}"

onsuccess:
	print("Workflow finished without errors")

onerror:
	print("An error occurred during runtime")

onstart:
	print("Setting up and running pipeline")
