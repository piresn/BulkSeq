# BulkSeq

<b>1.snpFile.R:</b> retrieve publicly available snp data for the Cvi-0 and Ler-1 accessions of <i>Arabidopsis thaliana</i>, merge and output a reformatted <i>snpm.csv</i> file

<b>2.cleanCounts.R:</b> merge information from <i>snpm.csv</i> file with an <i>out.vcf</i> file, filters and outputs a <i>counts.csv</i> file with allele counts

<b>3.pool.R:</b> combine allele frequencies from two samples (obtained with <i>cleanCounts.R</i>) and calculate relative frequencies along chromosomes

#Requires:
R version 3.3.1, scales_0.4.0, ggplot2 2.1.0, zoo 1.7-13
