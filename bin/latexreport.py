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
    results = {}
    for f in files:
        file_id = os.path.splitext(os.path.basename(f))[0]
        try:
            with open(f) as infile:
                line = infile.readline().strip()
                # Expect format like: "('stx2a_1_1', 0.0)"
                match = re.match(r"\('([^']+)',\s*([\d.]+)\)", line)
                if match:
                    allele, dist = match.groups()
                    results[file_id] = {"Allele": allele, "Distance": dist}
                else:
                    results[file_id] = {"Allele": "N/A", "Distance": "N/A"}
        except Exception as e:
            results[file_id] = {"Allele": f"Error: {e}", "Distance": "N/A"}
    return results


def parse_motif(files):
    results = {}
    for f in files:
        file_id = os.path.splitext(os.path.basename(f))[0]
        try:
            with open(f) as infile:
                last_line = infile.readlines()[-1].strip()
                match = re.search(r"subtype \['([a-z])'\]", last_line)
                if match:
                    results[file_id] = match.group(1)
                else:
                    results[file_id] = "N/A"
        except Exception as e:
            results[file_id] = f"Error: {e}"
    return results


def add_table(doc, headers, rows):
    with doc.create(Tabular("|" + "c|" * len(headers))) as table:
        table.add_hline()
        table.add_row(headers)
        table.add_hline()
        for row in rows:
            table.add_row(row)
        table.add_hline()


def safe_read_tsv(filepath, required_cols):
    """Read TSV and return only required columns, handling missing cols gracefully."""
    try:
        df = pd.read_csv(filepath, sep="\t")
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            # fill with 'N/A' if column missing
            for col in missing:
                df[col] = "N/A"
        return df[required_cols]
    except Exception as e:
        # Return dummy dataframe if file can't be read
        return pd.DataFrame([{col: f"Error: {e}" for col in required_cols}])


def main():
    parser = argparse.ArgumentParser(description="Generate PDF summary report with pylatex.")
    parser.add_argument("--closestleaf", nargs="*", help="Closestleaf text files")
    parser.add_argument("--motif", nargs="*", help="Motif text files")
    parser.add_argument("--stxtype", help="stx_type tab-delimited file")
    parser.add_argument("--blast", help="blast tab-delimited file")
    parser.add_argument("--output", default="report", help="Output PDF filename (without extension)")
    args = parser.parse_args()

    doc = Document()
    doc.preamble.append(Command("title", "Summary Report"))
    doc.preamble.append(Command("author", "Automated Pipeline"))
    doc.preamble.append(Command("date", NoEscape(r"\today")))
    doc.append(NoEscape(r"\maketitle"))

    # Section 1: Closestleaf + Motif
    if args.closestleaf or args.motif:
        with doc.create(Section("Closestleaf & Motif Results")):
            closest_data = parse_closestleaf(args.closestleaf) if args.closestleaf else {}
            motif_data = parse_motif(args.motif) if args.motif else {}

            # Collect union of all file IDs
            all_ids = sorted(set(closest_data.keys()) | set(motif_data.keys()))

            rows = []
            for fid in all_ids:
                allele = closest_data.get(fid, {}).get("Allele", "N/A")
                dist = closest_data.get(fid, {}).get("Distance", "N/A")
                motif_val = motif_data.get(fid, "N/A")
                rows.append([fid, allele, dist, motif_val])

            headers = ["File ID", "Allele", "Distance", "Motif"]
            add_table(doc, headers, rows)

    # Section 2: stx_type file
    if args.stxtype:
        with doc.create(Section("stx_type Results")):
            headers = ["stx_type", "operon", "identity"]
            df = safe_read_tsv(args.stxtype, headers)
            rows = df.values.tolist()
            add_table(doc, headers, rows)

    # Section 3: blast file
    if args.blast:
        with doc.create(Section("Blast Results")):
            headers = ["Template", "Template_Identity", "Template_Coverage"]
            df = safe_read_tsv(args.blast, headers)
            rows = df.values.tolist()
            add_table(doc, headers, rows)

    # Generate PDF
    doc.generate_pdf(args.output) #, clean_tex=False)
    # doc.generate_pdf(args.output, compiler="tectonic", command=["tectonic"], clean_tex=False)

if __name__ == "__main__":
    main()
