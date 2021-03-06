#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import pandas as pd
import shutil
import subprocess


parser = argparse.ArgumentParser()
parser.add_argument('--sequences', help='FASTA file of reference sequences')
parser.add_argument('--sampleName')
parser.add_argument('--assembly')
parser.add_argument('--minLength', type=int)
parser.add_argument('--default', help='reference sequence if no hits')
args = parser.parse_args()

sequence = SeqIO.read(args.assembly, 'fasta')
# check that the query has enough sequence to BLAST
if len(sequence)-sequence.seq.count('N') > args.minLength/2:
    subprocess.run(" ".join(["blastn", "-db", "blast_seqs.nt", "-query",
                    args.assembly, "-num_threads", "32", "-out",
                    f"{args.sampleName}.blast.tsv", "-outfmt", "'6",
                    "sacc", "nident", "pident", "length", "mismatch",
                    "gapopen", "qstart", "qend", "sstart", "send",
                    "evalue", "bitscore'"]), shell=True)

    df = pd.read_csv(f'{args.sampleName}.blast.tsv', sep='\t', names=['sacc', 'nident', 'pident',
                                                          'length', 'mismatch', 'gapopen', 'qstart',
                                                          'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    if df['length'].max() < args.minLength/4:
        # One HSP should extend to at least a few thousand bp
        # since we only expect a few SNPs
        # Otherwise the alignment is probably bad
        shutil.copyfile(args.default, f'{args.sampleName}_nearest_blast.fasta')
    else:
        try:
            # In case the TSV is empty
            top_hit = df[df['bitscore']==df['bitscore'].max()]['sacc']
            sequences = SeqIO.index(args.sequences, 'fasta')
            top_hit_seqs = [sequences[h] for h in top_hit]
            SeqIO.write(top_hit_seqs, f'{args.sampleName}_nearest_blast.fasta', 'fasta')
        except IndexError:
            shutil.copyfile(args.default, f'{args.sampleName}_nearest_blast.fasta')
else:
    shutil.copyfile(args.default, f'{args.sampleName}_nearest_blast.fasta')
    with open(f'{args.sampleName}.blast.tsv', 'w+') as f:
        f.write('Empty query sequence.')
