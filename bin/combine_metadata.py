#!/usr/bin/env python3

import argparse
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('--nextstrain_metadata', help='nextstrain metadata TSV')
parser.add_argument('--sample_metadata', help='sample metadata TSV, must have strain as a column')
args = parser.parse_args()

nextstrain_meta = pd.read_csv(args.nextstrain_metadata, sep='\t')
sample_meta = pd.read_csv(args.sample_metadata, sep='\t')

df = pd.concat([nextstrain_meta, sample_meta],
                axis=0,
                join='outer',
                sort=False)

# Drop duplicate strains and only keep the first row
df = df.drop_duplicates(subset='strain')
df = df.fillna('?')
df.to_csv('metadata.tsv', sep='\t', index=False)