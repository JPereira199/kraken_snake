#!/usr/bin/env bash
#set -euo pipefail

# Adjust these paths if needed:
RAW_DIR="/media/Disco11/pperaza/pperaza/pperaza-novogene/metagem_data/dataset"
OUT_DIR="../test_data/subsample_1k"
mkdir -p "$OUT_DIR"

# Number of reads to keep:
N=1000
# FASTQ: 4 lines per read
#LINES=$(( N * 4 ))

# Loop through each sample directory
for sample_name in $(ls "$RAW_DIR"); do
  sample_path="$RAW_DIR/$sample_name"
  [ -d "$sample_path" ] || continue

  mkdir -p "$OUT_DIR/$sample_name"

  # Find read files (_1 and _2 paired-end files)
  read1=$(find "$sample_path" -type f -regex ".*_R?1\.\(fq\|fastq\)\.gz")
  read2=$(find "$sample_path" -type f -regex ".*_R?2\.\(fq\|fastq\)\.gz")

  echo Read1: $read1

  if [[ -f "$read1" ]]; then
    outfile1="$OUT_DIR/$sample_name/${sample_name}_1.fq.gz"
    echo "SeqKit head $read1 → $outfile1"
    seqkit head -n "$N" -w 0 "$read1" -o "$outfile1"
  fi

  if [[ -f "$read2" ]]; then
    outfile2="$OUT_DIR/$sample_name/${sample_name}_2.fq.gz"
    echo "SeqKit head $read2 → $outfile2"
    seqkit head -n "$N" -w 0 "$read2" -o "$outfile2"
  fi
done
