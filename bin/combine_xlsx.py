#!/usr/bin/env python3
# Merges per-sample TSV rows (produced by xlsxreport.py) into one Excel file.

import argparse
import pandas as pd
import json

def parse_parameters(params_path):
    with open(params_path) as f:
        params_dict = json.load(f)
    manifest = params_dict.pop("manifest", {})
    df = pd.DataFrame.from_dict(params_dict, orient='index', columns=["value"])
    return df, manifest.get("homePage", ""), manifest.get("version", "")

def parse_versions(versions_path):
    """
    Parse versions.json (a flat list of [process, tool, version] triples,
    written by COMBINE_XLSX via topic channel + JsonBuilder) into a
    deduplicated, sorted DataFrame.
    """
    with open(versions_path) as f:
        versions_list = json.load(f)
    versions_dict = {it[0]: it[1:] for it in versions_list}
    return versions_dict


def infotab(toolversions, plink, pversion):
    stxtyperv = toolversions.pop("STXTYPER", ['','no_version'])[1]
    kmav = toolversions.pop("KMA", ['', 'no_version'])[1]
    rows = [f"This Nextflow analysis pipeline ({plink}, {pversion}) and database (https://github.com/OLC-Bioinformatics/STxDB) identifies and characterizes stx1/stx2 from assembly or read data, and outputs a results summary that includes:",
            "1. HC-Stx, which identifies full-length stx hits to an stxDB using assembled genomes. HC-Stx used the following external tools for its analysis:"]
    for key in toolversions:
        rows.append(f"\t{toolversions[key][0]} {toolversions[key][1]}")
    rows.extend([f"2. NCBI-Stx identifies full-length and partial stx hits to an stxDB using assembled genomes.  NCBI-stx  (https://github.com/ncbi/stxtyper {stxtyperv}).",
            f"3. KMA-Stx identifies full-length and partial stx hits to an stxDB using genomic reads. KMA-stx  (https://bitbucket.org/genomicepidemiology/kma/src/master {kmav}).",
            "",
            "If any of these analyses fail, no summary report is generated. Blank fields in a report indicate that there was no hit identified for that analysis."])
    df = pd.DataFrame(rows, columns=[''])
    return df


def makepretty(ws):
    for col_cells in ws.columns:
        max_len = max(len(str(c.value)) if c.value is not None else 0
                for c in col_cells)
        ws.column_dimensions[col_cells[0].column_letter].width = max_len + 2

        
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tsvs", nargs="+", help="Per-sample TSV files")
    parser.add_argument("--versions", help="YAML formatted tool versions")
    parser.add_argument("--params", help="YAML formatted pipeline parameters")
    parser.add_argument("--output", default="summary_report.xlsx")
    args = parser.parse_args()

    frames = []
    for f in args.tsvs:
        try:
            frames.append(pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False))
        except Exception as e:
            print(f"Warning: could not read {f}: {e}", flush=True)

    if not frames:
        raise SystemExit("No TSV files could be read — aborting.")

    df = pd.concat(frames, ignore_index=True).sort_values("Sequence ID")
    
    if args.versions:
        vdict = parse_versions(args.versions)
    else:
        vdict = {}

    if args.params:
        paramdf, plink, pversion = parse_parameters(args.params)
    else:
        plink = ""
        pversion = ""

    with pd.ExcelWriter(args.output, engine="openpyxl") as writer:
        infotab(vdict, plink, pversion).to_excel(writer, index=False, sheet_name="Info")
        wsi = writer.sheets["Info"]
        makepretty(wsi)

        df.to_excel(writer, index=False, sheet_name="Summary")
        # Auto-size columns for readability
        wss = writer.sheets["Summary"]
        makepretty(wss)

        if args.params:
            paramdf.to_excel(writer, index=True, sheet_name="Parameters")
            wsp = writer.sheets["Parameters"]
            makepretty(wsp)

    print(f"Wrote {args.output} ({len(df)} samples)", flush=True)


if __name__ == "__main__":
    main()
