#!/usr/bin/env python3
# Julie Shay and ChatGPT
# September 8, 2025
# Creates a PDF report from various pipeline output files
import argparse
import os
import re
import pandas as pd
from pylatex import Document, Section, Tabular, Command
from pylatex.utils import NoEscape


def parse_closestleaf(files):
    results = {"Allele": [], "Distance": []}
    for f in files:
        line = ""
        try:
            with open(f) as infile:
                line = infile.readline().strip()
        except:
            continue
        if line:
            # Expect format like: "('stx2a_1_1', 0.0)"
            match = re.match(r"\('([^']+)',\s*([\d.]+)\)", line)
            if match:
                allele, dist = match.groups()
                results["Allele"].append(allele)
                results["Distance"].append(dist)
            else:
                results["Allele"].append("")
                results["Distance"].append("")
        else:
            results["Allele"].append("File error")
            results["Distance"].append("File error")
    
    return results


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


def add_table(doc, headers, rows):
    with doc.create(Tabular("|" + "c|" * len(headers))) as table:
        table.add_hline()
        table.add_row(headers)
        table.add_hline()
        for row in rows:
            table.add_row(row)
        table.add_hline()


def main():
    parser = argparse.ArgumentParser(description="Generate PDF summary report with pylatex.")
    parser.add_argument("--id", type=str, required=True, help="seqID")
    parser.add_argument("--stx", type=str, nargs="*", help="stx type")
    parser.add_argument("--contig", type=str, nargs="*", help="contig name")
    parser.add_argument("--loc", type=str, nargs="*", help="location on contig")
    parser.add_argument("--closestleaf", type=str, nargs="*", help="Closestleaf text files")
    parser.add_argument("--motif", type=str, nargs="*", help="Motif text files")
    parser.add_argument("--stxtyper", type=str, help="stx_type tab-delimited file")
    parser.add_argument("--kma", type=str, help="kma tab-delimited file")
    parser.add_argument("--output", type=str, default="report", help="Output PDF filename (without extension)")
    args = parser.parse_args()

    doc = Document()
    doc.preamble.append(Command("title", "Stx Analysis for " + args.id))
    doc.preamble.append(Command("date", NoEscape(r"\today")))
    doc.append(NoEscape(r"\maketitle"))

    # Section 1: Closestleaf + Motif
    if args.closestleaf or args.motif:
        with doc.create(Section("Stx-subtyper Results")):
            closest_data = parse_closestleaf(args.closestleaf) if args.closestleaf else {}
            motif_data = parse_motif(args.motif, args.stx) if args.motif else {}
            headers = ["Gene", "Allele", "Distance", "Motif", "Contig", "Location"]
            df = pd.DataFrame([args.stx, closest_data["Allele"], closest_data["Distance"], motif_data, args.contig, args.loc])
            # Collect union of all file IDs
            rows = df.T.values.tolist()
            print(rows)
            add_table(doc, headers, rows)

    # Section 2: stx_type file
    if args.stxtyper:
        with doc.create(Section("Stxtyper Results")):
            cols = ["stx_type", "identity", "operon", "target_contig", "target_start", "target_stop"]
            try:
                df = pd.read_csv(args.stxtyper, sep="\t", usecols=cols, dtype=str)[cols]
                df["Location"] = df["target_start"] + "-" + df["target_stop"]
                df = df.drop(columns=["target_start", "target_stop"])
                headers = ["Stx Type", "Identity", "Operon", "Contig", "Location"]
                rows = df.values.tolist()
                add_table(doc, headers, rows)
            except:
                doc.append("There was a problem reading the stxtyper file.")

    # Section 3: KMA file
    if args.kma:
        with doc.create(Section("KMA Results")):
            try:
                cols = ["#Template", "Query_Identity", "Depth"]
                df = pd.read_csv(args.kma, sep="\t", usecols=cols, dtype=str)[cols]
                df["#Template"] = df["#Template"].str.replace(r"\|.*", "", regex=True)
                headers = ["Allele", "Identity", "Depth"]
                rows = df.values.tolist()
                add_table(doc, headers, rows)
            except:
                doc.append("There was a problem reading the kma output file.")

    # Generate PDF
    doc.generate_pdf(args.output) #, clean_tex=False)

if __name__ == "__main__":
    main()
