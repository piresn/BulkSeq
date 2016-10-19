#!/usr/bin/env python3
# -*- coding: utf-8 -*-

print("""

    *** Usage ./mask.py (r/m) FASTA_file SNP_file ***

mask.py returns a FASTA file with the SNPs indicated in SNP_file masked
SNP_file should be csv with: chromosome, position, ref_base, new_base

- m option: mask SNPs with Ns
- r option: replace bases as in SNP_file

""")

import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

################################################
######## import snp data
################################################

snpfile = open(sys.argv[3])
snp = csv.reader(snpfile)

info = []
for row in snp:
    info.append(row)

################################################
# Define functions
################################################

#subset snp file per chromsome
def subinfo(chr, info):

    subset = []

    for i in range(len(info)):
        if info[i][0]==chr:
            subset.append(info[i])
    return(subset)


# replace chromosome sequence (index 0) according to snp file (index (1)
def replace(seq, info, chr):

    flag=[]
    out = seq
    subset = subinfo(chr, info)

    for i in range(len(subset)):

        #check if reference matches between fasta and in snp file
        if seq[int(subset[i][1])-1] != subset[i][2]:
            flag.append(i)

        # replace pos i-1 (index 0) with new base
        if sys.argv[1] =='r':
            out[int(subset[i][1])-1] = subset[i][3]

        # replace pos i-1 (index 0) with N
        if sys.argv[1] =='n':
            out[int(subset[i][1])-1] = 'N'


    print('Replaced', len(subset), 'base(s)')

    if len(flag)>0:
        print('Warning: ', len(flag), ' base(s) do not match between fasta and SNP info file.')

    return(out)


# main function, takes record.seq, creates list, applies replace function and returns new record.seq object
def main(record, info):

    print("\nProcessing", record.id)

    repl = replace(list(record.seq), info, record.id)
    out = SeqRecord(Seq("".join(repl)), id=record.id, description="")

    return(out)

################################################
# run
################################################

f = open(sys.argv[2])

# iterate over fasta sequences, apply main function and returns a generator
out = (main(record, info) for record in SeqIO.parse(f, 'fasta'))

#save to file
output_handle = open("out.fas", "w")
SeqIO.write(out, output_handle, "fasta")
output_handle.close()
