#!/usr/bin/env python3
# Julie Shay
# February 18, 2025
# Given BLAST outputs for each subunit, this script will see whether there is a set of A and B subunits
# who are on the same contig, close enough to each other, and in the same orientation.
# if so, it will print a fasta of the concatenated sequence

# also, maybe if A and B are exact matches to the same allele, it can say so?

import argparse
import re
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('-a', required=True, help="""blast output for subunit A""")
parser.add_argument('-b', required=True, help="""blast output for subunit B""")
parser.add_argument('-c', dest='cutoff', type=int, default=1500, help="""cutoff for length of gap between subunits""")
parser.add_argument('-l', dest='minlen', type=float, default=0.99, help="""cutoff for length of alignment relative to gene length (max 1)""")
parser.add_argument('-o', dest='oprefix', default='test', help="""prefix for output files""")
args = parser.parse_args()

connecter = "NNNNNNNNNNNN"

# remove redundancy: if hits overlap, only keep the top one
def remove_overlapping(df):
    df.drop_duplicates(subset=['qacc', 'qstart', 'qend'], keep='first', inplace=True)
    kept_rows = []  # Store row indices to keep
    
    for idx, row in df.iterrows():
        qstart, qend, qacc = row['qstart'], row['qend'], row['qacc']
        
        # Check if this row overlaps with any previously kept row with the same subject
        overlap = any(
            (qstart <= df.loc[kept_idx, 'qend'] and qend >= df.loc[kept_idx, 'qstart'])  
            for kept_idx in kept_rows if df.loc[kept_idx, 'qacc'] == qacc  # Same subject
        )

        if not overlap:
            kept_rows.append(idx)  # Store only the index

    return df.loc[kept_rows]  # Keep only the selected rows


# for a blast match to subunit A, look for an appropriate pair from the matches to subunit B
# assuming genome is query and subunits are subject
def parseb(arow, wholeb):
    # filter for the same contig, same direction
    subb = wholeb[(wholeb['qacc'] == arow["qacc"]) & ((wholeb['sframe'] > 0) == (arow["sframe"] > 0))]

    if not subb.empty:
        if arow['sframe'] > 0:
            # B should follow A
            bmatch = subb[(subb['qstart'] > arow['qend']) & (subb['qstart'] <= (arow['qend'] + args.cutoff))]
        else:
            # A should follow B
            bmatch = subb[(subb['qend'] < arow["qstart"]) & (subb['qend'] >= (arow["qstart"] - args.cutoff))]
        if not bmatch.empty:
            # changing the column names so they don't conflict with each other
            bmatch.columns = [col + "_B" for col in bmatch.columns]
            acopy = arow.copy(deep=True)
            acopy.index = [str(idx) + "_A" for idx in acopy.index]
            # merging the result
            bmatch[acopy.index] = acopy.values
            return bmatch


# For a matching subunit pair, checks whether they are perfect matches to the A and B subunits of the same allele
# I'm not sure whether this even makes sense with the way the reference allele sequences are organized
'''
def checkperfect(combinedrow):
    if ((combinedrow['mismatch_A'] == 0) and (combinedrow['mismatch_B'] == 0)) and ((combinedrow['length_A'] == combinedrow['slen_A') and (combinedrow['length_B'] == combinedrow['slen_B'])):
        aallele = re.search('_(.*)\|', combinedrow['qacc_A']).group(1)
        ballele = re.search('_(.*)\|', combinedrow['qacc_B']).group(1)
        if aallele == ballele:
            return aallele
        else:
            return ""
    else:
        return ""

'''

def concatenateseqs(combinedrow):
    if combinedrow['sframe_A'] > 0:
        for i in range(0, ((combinedrow['sstart_A'] - 1) % 3)):
            combinedrow['qseq_A'] = 'N' + combinedrow['qseq_A']
        while (len(combinedrow['qseq_A']) % 3):
            combinedrow['qseq_A'] += 'N'
        for i in range(0, ((combinedrow['sstart_B'] - 1) % 3)):
            combinedrow['qseq_B'] = 'N' + combinedrow['qseq_B']
        while (len(combinedrow['qseq_B']) % 3):
            combinedrow['qseq_B'] += 'N'
        combinedseq = Seq(combinedrow['qseq_A'] + connecter + combinedrow['qseq_B'])
        description = "-".join((str(combinedrow['qstart_A']), str(combinedrow['qend_B'])))
    else:
        for i in range(0, ((combinedrow['send_A'] - 1) % 3)):
            combinedrow['qseq_A'] += 'N'
        while (len(combinedrow['qseq_A']) % 3):
            combinedrow['qseq_A'] = 'N' + combinedrow['qseq_A']
        for i in range(0, ((combinedrow['send_B'] - 1) % 3)):
            combinedrow['qseq_B'] += 'N'
        while (len(combinedrow['qseq_B']) % 3):
            combinedrow['qseq_B'] = 'N' + combinedrow['qseq_B'] 
        combinedseq = Seq(combinedrow['qseq_B'] + connecter + combinedrow['qseq_A']).reverse_complement()
        description = "-".join((str(combinedrow['qstart_B']), str(combinedrow['qend_A'])))
    
    return SeqRecord(combinedseq.replace("-", ""), id=combinedrow['qacc_A'], description=description)


adf = pd.read_csv(args.a, sep="\t")
bdf = pd.read_csv(args.b, sep="\t")


# filter for alignments that are close to full length
adf.drop(adf[adf['length'] <= (adf['slen'] * args.minlen)].index, inplace=True)
bdf.drop(bdf[bdf['length'] <= (bdf['slen'] * args.minlen)].index, inplace=True)

# remove overlapping hits
adf = remove_overlapping(adf)
bdf = remove_overlapping(bdf)

# print these
adf.to_csv(args.oprefix + "_filtered_A.tsv", sep="\t", index=False)
bdf.to_csv(args.oprefix + "_filtered_B.tsv", sep="\t", index=False)


if adf.empty or bdf.empty:
    exit()

# look for appropriate a and b subunit pairs
results = adf.apply(parseb, args=(bdf,), axis=1).dropna().tolist()

if not results:
    exit()

matches = pd.concat(results, ignore_index=True)
# matches['allele'] = matches.apply(checkperfect, axis=1)

# print the results?
matches.to_csv(args.oprefix + "_properpairs.tsv", sep="\t", index=False)

# Produce a fasta file with the concatenated sequences for input to other scripts
nucseqs = matches.apply(concatenateseqs, axis=1)

SeqIO.write(nucseqs, args.oprefix + "_nuc.fasta", "fasta")
# Translate and collect protein sequences
protseqs = [nucseq.translate(id=nucseq.id, description=nucseq.description) for nucseq in nucseqs]

# Write all protein sequences at once
SeqIO.write(protseqs, args.oprefix + "_prot.fasta", "fasta")
