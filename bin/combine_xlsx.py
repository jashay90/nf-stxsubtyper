#!/usr/bin/env python3
# Merges per-sample TSV rows (produced by xlsxreport.py) into one Excel file.

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tsvs", nargs="+", help="Per-sample TSV files")
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

    with pd.ExcelWriter(args.output, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="Summary")

        # Auto-size columns for readability
        ws = writer.sheets["Summary"]
        for col_cells in ws.columns:
            max_len = max(len(str(c.value)) if c.value is not None else 0
                         for c in col_cells)
            ws.column_dimensions[col_cells[0].column_letter].width = min(max_len + 2, 60)

    print(f"Wrote {args.output} ({len(df)} samples)", flush=True)


if __name__ == "__main__":
    main()
