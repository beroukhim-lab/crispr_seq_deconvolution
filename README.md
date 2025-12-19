# crispr_seq_deconvolution
This repo describes an analysis pipeline for local analysis of CRISPR-Seq data using crispresso.


### Step 1: Split input fastq file
In my wetlab workflow I typically pool multiple PCR reactions (with different primer pairs) then sequence them in the same run. We want to deconvolute these and create individual fastq files for each primer pair. The first step in this analysis pipeline is to take a single fastq file and split it using a list of primer pairs.
```
split_fastq.py --fastq input.fastq --annotation annotation.txt --outdir /path/to/output/dir/
```
--annotation uses a standard input file used in multiple steps in this pipeline. Not all fields are directly used in split_fastq.py. File structure:
Column 1: User defined primer set ID. This id will be appended to output files
Column 2: First primer sequence (primer order does not matter, all orientations will be checked)
Column 3: Second primer sequence (primer order does not matter, all orientations will be checked)
Column 4: Guide sequence (if using Cas12a, each guide should be on a different line)
Column 5: Complete amplicon sequence, including primer sequences
It is assumed that the annotation file is in a tab delimited format.


### Step 2: Construct a read count matrix to check for large cross contamination or PCR failures
```
qual_check.py --annotation annotation.txt --input_dir
```
--annotation is the same annotation file as above
--input_dir is a path to the directory where the fastq files are stored. fastq files which do not begin with a primer set ID in the annotation file will be ignored. File name must start with 


This script will export a read count df and a read count matrix. You should manually check the read count matrix to ensure that there are no sample swaps. If you have a large number of samples and there is a pattern to your primer and fastq names, this may be faster by writing a custom script to analyze the output df.


### Step 3: Run CRISPresso on all output files
```
run_crispresso.sh --annotation annotation.txt --input_dir
```


