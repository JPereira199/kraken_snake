#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path
import argparse
from itertools import chain

import subprocess
import sys

def run_command(cmd, stop_error_str: bool = True):
    """Execute a command and check for stderr containing 'error' or 'ERROR'."""
    print(f" Running: {' '.join(cmd)}")
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,  # ensures output is returned as strings, not bytes
        check=False  # we will handle errors manually
    )

    # Print captured outputs for debugging/logging
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)

    if stop_error_str:
        # Check for explicit error in stderr
        if "error" in result.stderr.lower():
            print("❌ Bracken command failed: found 'error' in stderr output.")
            sys.exit(1)
        if "error" in result.stdout.lower():
            print("❌ Bracken command failed: found 'error' in stdout output.")
            sys.exit(1)


    # Optionally also check return code (if needed)
    if result.returncode != 0:
        print(f"❌ Command exited with non-zero code: {result.returncode}")
        sys.exit(result.returncode)

    return result.returncode


# ---------------------------- CLI ARGUMENTS ----------------------------
parser = argparse.ArgumentParser(
    description="Estimate species abundances using Bracken from Kraken2 report files"
)
parser.add_argument(
    "--input-kraken-dir",
    required=True,
    type=Path,
    help="Parent folder containing one subdirectory per sample (each with a *.kraken2.report file)"
)
parser.add_argument(
    "--bracken-db",
    required=True,
    type=Path,
    help="Path to the Bracken database folder (will build if not already built)"
)
parser.add_argument(
    "--read-length",
    type=int,
    default=100,
    help="Read length used for Bracken estimation (default: 100)"
)
parser.add_argument(
    "--kraken-db-kmer",
    type=int,
    default=35,
    help="Kmer size used to build the Kraken2 database"
)
parser.add_argument(
    "--level",
    choices=["D", "P", "C", "O", "F", "G", "S"],
    default="S",
    help="Taxonomic level for Bracken output (D=Domain, P=Phylum, C=Class, O=Order, F=Family, G=Genus, S=Species; default: S)"
)
parser.add_argument(
    "--threshold",
    type=int,
    default=0,
    help="Number of reads requiered PRIOR to estimate the abundace (default: 0)"
)
parser.add_argument(
    "--threads",
    type=int,
    default=10,
    help="Number of threads used to build bracken database"
)

parser.add_argument(
    "--output-dir",
    required=True,
    type=Path,
    help="Directory to write Bracken results (one subfolder per sample)"
)
args = parser.parse_args()

input_kraken_dir = args.input_kraken_dir
bracken_db = args.bracken_db
read_length = args.read_length
kraken_db_kmer = args.kraken_db_kmer
tax_level = args.level
threshold = args.threshold
threads = args.threads
output_dir = args.output_dir

# Validate input directories
if not input_kraken_dir.exists():
    print(f"⚠️  Input Kraken directory not found: {input_kraken_dir}")
    raise SystemExit(1)
# Ensure Bracken DB directory exists (or create it if building from scratch)
bracken_db.mkdir(parents=True, exist_ok=True)

# ---------------------------- BUILD Bracken DB IF NEEDED ----------------------------
# Check for any existing k-mer distribution files (*.kmer_distrib*)
existing_kmer = list(bracken_db.glob("*.kmer_distrib*"))
if not existing_kmer:
    print(f"📦 Building Bracken database in: {bracken_db}")
    # Example k-mer length of 35; adjust if needed
    kmer_length = kraken_db_kmer
    bracken_build_cmd = [
        "bracken-build",
        "-d", str(bracken_db),
        "-t", str(threads),
        "-k", str(kmer_length),
        "-l", str(read_length)
    ]
    
    try:
        run_command(bracken_build_cmd)
        print(f"✅ Bracken database built successfully in {bracken_db}")
    except Exception as e:
        print(f"❌ Bracken-build raised: {type(e).__name__}: {e}")
        raise SystemExit(1)

else:
    print(f"✅ Found existing Bracken k-mer distribution files in {bracken_db}, skipping build.")

# Ensure output directory exists
output_dir.mkdir(parents=True, exist_ok=True)

# ----------------------------- MAIN LOOP -----------------------------
failed_samples = []

for sample_dir in sorted(input_kraken_dir.iterdir()):
    if not sample_dir.is_dir():
        continue
    sample_name = sample_dir.name
    print(f"\n🔬 Processing sample: {sample_name}")

    # Find the Kraken2 report file: *.kraken2.report
    kraken_reports = list(sample_dir.glob("*.kraken2.report"))
    if len(kraken_reports) != 1:
        print(f"⚠️  Could not find exactly one Kraken2 report (*.kraken2.report) in {sample_dir}")
        failed_samples.append(sample_name)
        continue

    report_file = kraken_reports[0]

    # Prepare output directory for this sample
    sample_out = output_dir / sample_name
    sample_out.mkdir(parents=True, exist_ok=True)

    # Define Bracken output files
    bracken_out = sample_out / f"{sample_name}.bracken.{tax_level}.txt"
    bracken_cov = sample_out / f"{sample_name}.bracken.{tax_level}.cov.txt"

    # Build Bracken command (abundance estimation)
    bracken_cmd = [
        "bracken",
        "-d", str(bracken_db),
        "-i", str(report_file),
        "-o", str(bracken_out),
        "-r", str(read_length),
        "-l", str(tax_level),
        "-t", str(threshold)
    ]

    # Execute Bracken abundance estimation
    try:
        run_command(bracken_cmd)
        print(f"✅ Bracken abundance estimation complete for {sample_name}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Bracken failed for {sample_name}: {e}")
        failed_samples.append(sample_name)
        continue

# Final status report
if failed_samples:
    print("\n❌ Bracken estimation failed for the following samples:")
    for sample in failed_samples:
        print(f" - {sample}")
    raise SystemExit(1)
else:
    print("\n✨ Bracken estimation completed successfully for all samples!")
