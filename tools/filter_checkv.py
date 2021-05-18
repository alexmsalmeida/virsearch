#!/usr/bin/env python

import sys
import Bio
from Bio import SeqIO
from Bio import SeqRecord

if len(sys.argv) < 3:
    print("usage: script.py in.fa checkv_summary.tsv out.fa")
    sys.exit(1)

contigs = set()
with open(sys.argv[2]) as f:
    for line in f:
        line = line.rstrip()
        cols = line.split("\t")
        if cols[0] != "contig_id":
            contig_id = cols[0].split("_length")[0]
            contigs.add(contig_id)

with open(sys.argv[1]) as f, open(sys.argv[3], "w") as fout:
    for record in SeqIO.parse(f, "fasta"):
        contig_id = record.id.split("_length")[0]
        if contig_id in contigs and len(record.seq) >= 10000:
            SeqIO.write(record, fout, "fasta")
            
