#!/bin/bash

#Check read quality with FastQC
fastqc input.fastq

# Remove residual adaptors, trim and/or discard low quality reads using cutadapt (Martin 2011).
# E.g trim low quality 3’ ends and remove contaminating Ilumina TruSeq adaptors:
cutadapt --minimum-length 20 -q 10 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o trim.fastq input.fastq


# Recheck read quality using FastQC and modify previous step if required


#	Retrieve a reference genome sequence (e.g. the TAIR10 genome release from www.arabidopsis.org) and index using Bowtie2
bowtie2-build -f TAIR10.fa index

# Align reads to the indexed genome using Bowtie2 and the SAMtools package
bowtie2 -x index -U trim.fastq | samtools view -b - > mapped.bam

# Index and sort alignment, and call SNPs with mpileup and bcftools
samtools sort mapped.bam mappednsorted
samtools index mappednsorted.bam
samtools mpileup -uf TAIR10.fa mappednsorted.bam | bcftools call -mv - > out.vcf
