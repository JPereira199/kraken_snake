#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path
import argparse
from itertools import chain

def run_command(cmd):
    """Execute a command and check for success."""
    print(f" Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=True)
    return result.returncode

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected (True/False).")

# ---------------------------- CLI ARGUMENTS ----------------------------
parser = argparse.ArgumentParser(
    description="Run Kraken2 on a set of paired metagenomic samples")
parser.add_argument("--input-raw-dir",required=True, type=Path, help="Parent folder containing one subdirectory per sample (each with _1 and _2 FASTQ files)")
parser.add_argument("--params-kraken2-db",required=True,type=Path,help="Path to Kraken2 database")
parser.add_argument("--params-memory-mapping",type=str2bool,default=False,help="Use 'True' to avoid loading DB into RAM (slower but uses less memory).")
parser.add_argument("--params-threads",type=int,default=8,help="Number of threads to use for Kraken2 (default: 8)")
parser.add_argument("--output-dir",required=True,type=Path,help="Directory to write Kraken2 results (one subfolder per sample)")
args = parser.parse_args()

input_raw_dir = Path(args.input_raw_dir)
kraken2_db = Path(args.params_kraken2_db)
memory_mapping = args.params_memory_mapping
cpus = args.params_threads
output_dir = Path(args.output_dir)

# Ensure the Kraken2 database path exists
if not kraken2_db.exists():
    print(f"⚠️  Kraken2 database not found: {kraken2_db}")
    raise SystemExit(1)

# Ensure output directory exists
output_dir.mkdir(parents=True, exist_ok=True)

# ----------------------------- MAIN LOOP -----------------------------
failed_samples = []

for sample_dir in sorted(input_raw_dir.iterdir()):
    if not sample_dir.is_dir():
        continue
    sample_name = sample_dir.name
    print(f"\n🔬 Processing sample: {sample_name}")

    # Find R1 and R2 reads by globbing separately for fq.gz and fastq.gz
    r1_list = list(chain(
        sample_dir.glob("*1.fq.gz"),
        sample_dir.glob("*1.fastq.gz"),
        sample_dir.glob("*1.trimmed.fq.gz"),
        sample_dir.glob("*1.subsampled.fq.gz")
    ))
    r2_list = list(chain(
        sample_dir.glob("*2.fq.gz"),
        sample_dir.glob("*2.fastq.gz"),
        sample_dir.glob("*2.trimmed.fq.gz"),
        sample_dir.glob("*2.subsampled.fq.gz")
    ))

    if len(r1_list) != 1 or len(r2_list) != 1:
        print(f"⚠️  Could not find exactly one pair of FASTQ files in {sample_dir}")
        failed_samples.append(sample_name)
        continue

    r1 = r1_list[0]
    r2 = r2_list[0]

    # Prepare output directory for this sample
    sample_out = output_dir / sample_name
    sample_out.mkdir(parents=True, exist_ok=True)

    # Paths for Kraken2 report and output
    report_file = sample_out / f"{sample_name}.kraken2.report"
    output_file = sample_out / f"{sample_name}.kraken2.output"

    # Build Kraken2 command
    kraken2_cmd = [
        "kraken2",
        "--db", str(kraken2_db),
        "--paired", str(r1), str(r2),
        "--threads", str(cpus),
        "--report", str(report_file),
        "--output", str(output_file)
    ]

    if memory_mapping:
        kraken2_cmd.append("--memory-mapping")

    # Execute Kraken2
    try:
        run_command(kraken2_cmd)
        print(f"✅ Kraken2 classification complete for {sample_name}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Kraken2 failed for {sample_name}: {str(e)}")
        failed_samples.append(sample_name)

# Final status report
if failed_samples:
    print("\n❌ Failed samples:")
    for sample in failed_samples:
        print(f" - {sample}")
    exit(1)
else:
    print("\n✨ All samples processed successfully!")
