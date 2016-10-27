# BulkSeq

This set of scripts is described in a chapter entitled 'Identification of parent-of-origin-dependent QTLs using bulk-segregant sequencing (Bulk-Seq)' that will appear in the Springer Protocol Series 'Methods in Molecular Biology' on Plant Chromatin in 2017.


<b>GenomeSNPmask.py: </b> Remove or replace known SNP positions from a genome sequence file (fasta)

<b>mapping.sh: </b> minimal set of commands to filter and map reads from a fastq file, call SNPs and output a <i>out.vcf</i> file with allele frequencies

<b>snpFile.R:</b> retrieve publicly available snp data for the Cvi-0 and Ler-1 accessions of <i>Arabidopsis thaliana</i>, merge and output a reformatted snp matrix (<i>snpm.txt</i>)

<b>cleanCounts.R:</b> merge information from the snp matrix (<i>snpm.txt</i>) with the measured allele frequencies (<i>out.vcf</i>), filters and outputs a <i>counts.csv</i> file with allele counts

<b>pool.R:</b> combine allele frequencies from two samples (obtained with <i>cleanCounts.R</i>) and calculate relative frequencies along chromosomes

<i>Requires:</i>

FastQC 0.11.3;
cutadapt 1.8.3;
Samtools 1.2 (using htslib 1.2.1);
Bowtie 2 2.2.9;
R 3.3.1;
scales_0.4.0;
ggplot2 2.1.0;
zoo 1.7-13;
Python 3.4.0;
Bio 1.65;

<HR>
<i> Example fastq datasets that can be used in this analysis are available in the ArrayExpress database (www.ebi.ac.uk/arrayexpress) under accession number E-MTAB-5196:

WT_pool_1 (1.56GB)

mea_pool_1 (2GB)</i>
