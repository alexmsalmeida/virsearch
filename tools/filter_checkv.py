#!/usr/bin/env python

import sys
import Bio
from Bio import SeqIO
from Bio import SeqRecord

if len(sys.argv) < 3:
    print("usage: script.py in.fa checkv_summary.tsv out.fa pro")
    sys.exit(1)

contigs = set()
with open(sys.argv[2]) as f:
    for line in f:
        line = line.rstrip()
        cols = line.split("\t")
        if cols[0] != "contig_id":
            if cols[9] != "NA":
                if int(cols[5]) > int(cols[6]) and int(cols[1]) >= 10000 and float(cols[9]) >= 50 and float(cols[11]) < 5 and float(cols[12]) <= 1.0 and "contig >1.5x" not in line:
                    contigs.add(cols[0])

with open(sys.argv[1]) as f, open(sys.argv[3], "w") as fout:
    for record in SeqIO.parse(f, "fasta"):
        if sys.argv[4] == "pro":
            contig_id = "_".join(record.id.split("_")[:-1])
        elif sys.argv[4] == "vir":
            contig_id = record.id.split()[0]
        if contig_id in contigs and len(record.seq) >= 10000:
            SeqIO.write(record, fout, "fasta")
            
