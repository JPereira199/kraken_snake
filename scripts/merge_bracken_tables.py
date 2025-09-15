#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import argparse
import sys

# ---------------------------- CLI ARGUMENTS ----------------------------
parser = argparse.ArgumentParser(
    description="Merge Bracken output tables from multiple samples into a single CSV per metric."
)
parser.add_argument(
    "--input-dir",
    required=True,
    type=Path,
    help="Directory containing one subdirectory per sample (each with a Bracken output file)"
)
parser.add_argument(
    "--level",
    choices=["D", "P", "C", "O", "F", "G", "S"],
    default="S",
    help="Taxonomic level used in Bracken output (D=Domain, ..., S=Species; default: S)"
)
parser.add_argument(
    "--metrics",
    nargs="+",
    default=["fraction_total_reads"],
    help="List of abundance metrics to include (e.g., fraction_total_reads new_est_reads)"
)
parser.add_argument(
    "--output-dir",
    type=Path,
    required=True,
    help="Directory to write the merged CSV tables"
)
args = parser.parse_args()

input_dir = args.input_dir
level = args.level
metrics_list = args.metrics
output_dir = args.output_dir

output_dir.mkdir(parents=True, exist_ok=True)

# -------------------------- Merge Bracken Tables --------------------------
bracken_dfs = {metric: [] for metric in metrics_list}

for sample_dir in sorted(input_dir.iterdir()):
    if not sample_dir.is_dir():
        continue

    sample_name = sample_dir.name
    bracken_file = sample_dir / f"{sample_name}.bracken.{level}.txt"

    if not bracken_file.exists():
        print(f"⚠️ Skipping {sample_name}: No Bracken file at {bracken_file}")
        continue

    try:
        df = pd.read_csv(bracken_file, sep="\t", header=0)
    except Exception as e:
        print(f"❌ Failed to read {bracken_file}: {e}")
        continue

    if "name" not in df.columns:
        print(f"⚠️ File {bracken_file} is missing required 'name' column.")
        continue

    df.set_index("name", inplace=True)

    for metric in metrics_list:
        if metric not in df.columns:
            print(f"⚠️ Metric '{metric}' not found in {bracken_file}, skipping.")
            continue

        metric_df = df[[metric]].rename(columns={metric: sample_name})
        bracken_dfs[metric].append(metric_df)

# Check if anything was loaded
if not any(bracken_dfs.values()):
    print("❌ No valid Bracken tables were loaded.")
    sys.exit(1)

# Merge and write each metric's table
for metric, df_list in bracken_dfs.items():
    if not df_list:
        print(f"⚠️ No data found for metric '{metric}', skipping export.")
        continue

    merged_df = pd.concat(df_list, axis=1).fillna(0)
    output_file = output_dir / f"bracken_merged_{level}_{metric}.csv"
    merged_df.to_csv(output_file)
    print(f"✅ Merged Bracken table saved to: {output_file}")
