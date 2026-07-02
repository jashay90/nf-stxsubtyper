#!/usr/bin/env python3
# Replaces latexreport.py
# Produces a single-row TSV (to stdout) summarising one sample.
# A separate combine_xlsx.py step assembles all rows into one Excel file.

import argparse
import re
import sys
import pandas as pd


SEP = ";"


def parse_closestleaf(files):
    alleles, dists = [], []
    for f in files:
        try:
            with open(f) as fh:
                line = fh.readline().strip()
        except Exception:
            alleles.append(""); dists.append("")
            continue
        m = re.match(r"\('([^']+)',\s*([\d.eE+\-]+)\)", line)
        if m:
            alleles.append(m.group(1))
            dists.append(m.group(2))
        else:
            alleles.append(""); dists.append("")
    return alleles, dists


def parse_motif(files, stxes):
    results = []
    twocount = 0
    for i in range(len(stxes)):
        if stxes[i] == "2":
            last_line = ""
            try:
                with open(files[twocount]) as infile:
                    last_line = infile.readlines()[-1].strip()
            except:
                continue
            if last_line:
                match = re.search(r"subtype \['([a-z])'\]", last_line)
                if match:
                    results.append(match.group(1))
                else:
                    results.append("")
            else:
                results.append("File error")
            twocount += 1
        else:
            results.append("N/A")
    return results

def join(lst):
    """Return semicolon-joined string, or empty string for None/empty."""
    if not lst:
        return ""
    return SEP.join(str(x) for x in lst if x is not None)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--id",          required=True)
    parser.add_argument("--stx",         nargs="*", default=[])
    parser.add_argument("--contig",      nargs="*", default=[])
    parser.add_argument("--loc",         nargs="*", default=[])
    parser.add_argument("--closestleaf", nargs="*", default=[])
    parser.add_argument("--motif",       nargs="*", default=[])
    parser.add_argument("--stxtyper",    default=None)
    parser.add_argument("--kma",         default=None)
    parser.add_argument("--samplename",  default=None)
    parser.add_argument("--coverage",    default=None)
    parser.add_argument("--output",      default="-",
                        help="Output TSV path, or '-' for stdout")
    args = parser.parse_args()

    row = {"Sequence ID": args.id}

    # sample name
    if args.samplename:
        row["Sample Name"] = args.samplename.replace('[', '').replace(']', '')

    # custom pipeline
    clean_stx = [s for s in (args.stx or []) if s not in ("[]", "null", "")]
    if clean_stx:
        alleles, dists = parse_closestleaf(args.closestleaf) if args.closestleaf else ([], [])
        motifs = parse_motif(args.motif, args.stx) if args.motif else ["N/A"] * len(args.stx)

        row["HC-Protein Allele"]    = join(alleles)
        row["HC-Distance"]  = join(dists)
        row["HC-Motif"]     = join(motifs)
        row["HC-Contig"]    = join(args.contig)
        row["HC-Location"]  = join(args.loc)
    else:
        row["HC-Protein Allele"] = row["HC-Distance"] = ""
        row["HC-Motif"]    = row["HC-Contig"] = row["HC-Location"] = ""

    # ncbi stxtyper
    if args.stxtyper:
        try:
            cols = ["stx_type", "identity", "operon",
                    "target_contig", "target_start", "target_stop"]
            df = pd.read_csv(args.stxtyper, sep="\t", usecols=cols, dtype=str)[cols]
            df["loc"] = df["target_start"] + "-" + df["target_stop"]
            row["NCBI-Stx Type"]     = join(df["stx_type"].tolist())
            row["NCBI-Identity"] = join(df["identity"].tolist())
            row["NCBI-Operon"]   = join(df["operon"].tolist())
            row["NCBI-Contig"]   = join(df["target_contig"].tolist())
            row["NCBI-Location"] = join(df["loc"].tolist())
        except Exception as e:
            for col in ("NCBI-Stx Type","NCBI-Identity","NCBI-Operon",
                        "NCBI-Contig","NCBI-Location"):
                row[col] = f"parse error: {e}"
    else:
        for col in ("NCBI-Stx Type","NCBI-Identity","NCBI-Operon",
                    "NCBI-Contig","NCBI-Location", "HC-Protein Allele",
                    "HC-Distance", "HC-Motif", "HC-Contig", "HC-Location"): 
            row[col] = "not_performed"

    # kma
    if args.kma:
        try:
            cols = ["#Template", "Query_Identity", "Depth"]
            df = pd.read_csv(args.kma, sep="\t", usecols=cols, dtype=str)[cols]
            df["#Template"] = df["#Template"].str.replace(r"\|.*", "", regex=True)
            row["KMA-Nucleotide Allele"]    = join(df["#Template"].tolist()).replace(" ", "")
            row["KMA-Identity"]  = join(df["Query_Identity"].tolist()).replace(" ", "")
            row["KMA-Depth"]     = join(df["Depth"].tolist()).replace(" ", "")
        except Exception as e:
            row["KMA-Nucleotide Allele"] = row["KMA-Identity"] = row["KMA-Depth"] = f"parse error: {e}"
    else:
        row["KMA-Nucleotide Allele"] = row["KMA-Identity"] = row["KMA-Depth"] = "not_performed"

    # depth of coverage
    if args.coverage:
        row["Assembly-Mean Coverage"] = args.coverage.replace('[', '').replace(']', '')

    # write row
    # Why doesn't fillna do anything?
    out_df = pd.DataFrame([row]).replace("", "no_result")
    if args.output == "-":
        out_df.to_csv(sys.stdout, sep="\t", index=False)
    else:
        out_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
