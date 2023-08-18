#!/usr/bin/env python

import argparse
import sys

def ren_fasta(args):
    for line in open(args.fasta_file, "r"):
        if line[0] == ">":
            name = line.strip("\n").replace(">","")
            print(">%s_%s" % (args.prefix, name))
        else:
            print(line.strip("\n"))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rename multifasta file')
    parser.add_argument('-f', dest='fasta_file', help='Input FASTA file')
    parser.add_argument('-p', dest='prefix', help='Header prefix')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()
        ren_fasta(args)
