#!/usr/bin/env python3
import subprocess
from pathlib import Path
import argparse
from itertools import chain
from concurrent.futures import ProcessPoolExecutor, as_completed
import tempfile
import shutil


def run_command(cmd):
    """Execute a command and check for success."""
    print(f" Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=True)
    return result.returncode


def process_sample(sample_dir, output_dir, subsample, threads):
    """
    Trim with fastp and then subsample paired FASTQ files using BBMap's reformat.sh
    via temporary files, writing exactly `subsample` read-pairs to output_dir/sample_name.
    Each job may use up to `threads` CPU cores.
    """
    sample_name = sample_dir.name
    print(f"\n🔬 Starting processing for sample: {sample_name}")

    # Locate R1/R2
    r1_list = list(chain(
        sample_dir.glob("*1.fq.gz"), sample_dir.glob("*1.fastq.gz")
    ))
    r2_list = list(chain(
        sample_dir.glob("*2.fq.gz"), sample_dir.glob("*2.fastq.gz")
    ))
    if len(r1_list) != 1 or len(r2_list) != 1:
        print(f"⚠️  Could not find a single pair of FASTQ files in {sample_dir}")
        return (sample_name, False)

    r1, r2 = r1_list[0], r2_list[0]

    # Prepare output folder
    sample_out = output_dir / sample_name
    sample_out.mkdir(parents=True, exist_ok=True)

    # Final subsampled outputs
    out_r1 = sample_out / f"{sample_name}_1.subsampled.fq.gz"
    out_r2 = sample_out / f"{sample_name}_2.subsampled.fq.gz"

    # Use a temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        try:
            # Copy raw reads to temp
            tmp_r1_raw = tmpdir / r1.name
            tmp_r2_raw = tmpdir / r2.name
            shutil.copy(r1, tmp_r1_raw)
            shutil.copy(r2, tmp_r2_raw)

            # Define trimmed outputs
            trim_r1 = tmpdir / f"{sample_name}_1.trimmed.fq.gz"
            trim_r2 = tmpdir / f"{sample_name}_2.trimmed.fq.gz"

            # Fastp log paths
            fastp_html = sample_out / f"{sample_name}.fastp.html"
            fastp_json = sample_out / f"{sample_name}.fastp.json"

            # Run fastp (uses up to `threads` cores)
            fastp_cmd = [
                "fastp",
                "-i", str(tmp_r1_raw),
                "-I", str(tmp_r2_raw),
                "-o", str(trim_r1),
                "-O", str(trim_r2),
                "-w", str(threads),
                "-h", str(fastp_html),
                "-j", str(fastp_json)
            ]
            run_command(fastp_cmd)
            print(f"✅ fastp completed for {sample_name}")

            # Define temp subsampled outputs
            tmp_out_r1 = tmpdir / f"{sample_name}_1.subsampled.tmp.fq.gz"
            tmp_out_r2 = tmpdir / f"{sample_name}_2.subsampled.tmp.fq.gz"

            # Run BBMap reformat.sh on trimmed reads (also uses up to `threads` cores)
            reformat_cmd = [
                "reformat.sh",
                f"in1={trim_r1}", f"in2={trim_r2}",
                f"out1={tmp_out_r1}", f"out2={tmp_out_r2}",
                f"reads={subsample}",
                f"threads={threads}"
            ]
            run_command(reformat_cmd)

            # Move final outputs into place
            shutil.move(str(tmp_out_r1), str(out_r1))
            shutil.move(str(tmp_out_r2), str(out_r2))
            print(f"✅ Subsampling completed for {sample_name}, outputs at {sample_out}")
            return (sample_name, True)
        except subprocess.CalledProcessError as e:
            print(f"❌ Pipeline failed for {sample_name}: {e}")
            return (sample_name, False)


def main():
    parser = argparse.ArgumentParser(
        description="Trim and subsample paired-end FASTQ files using fastp and BBMap's reformat.sh"
    )
    parser.add_argument(
        "--input-raw-dir", required=True, type=Path,
        help="Directory with one subdirectory per sample (each with *_1 and *_2 FASTQ files)"
    )
    parser.add_argument(
        "--params-subsample", required=True, type=int,
        help="Number of read-pairs to keep per sample"
    )
    parser.add_argument(
        "--output-dir", required=True, type=Path,
        help="Directory to write outputs (one subfolder per sample)"
    )
    parser.add_argument(
        "--threads", type=int, default=16,
        help="Cores per job (fastp and reformat.sh will use up to this many threads, default: 16)"
    )
    parser.add_argument(
        "--jobs", type=int, default=None,
        help="Number of samples to process in parallel (default: one job per sample)"
    )
    args = parser.parse_args()

    if not args.input_raw_dir.exists():
        print(f"⚠️  Input directory not found: {args.input_raw_dir}")
        raise SystemExit(1)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sample_dirs = sorted(d for d in args.input_raw_dir.iterdir() if d.is_dir())
    if not sample_dirs:
        print(f"⚠️  No sample subdirectories in: {args.input_raw_dir}")
        raise SystemExit(1)

    # If --jobs not set, run one job per sample
    jobs = args.jobs if args.jobs is not None else len(sample_dirs)

    failed = []
    with ProcessPoolExecutor(max_workers=jobs) as executor:
        futures = {
            executor.submit(
                process_sample,
                sample_dir,
                args.output_dir,
                args.params_subsample,
                args.threads
            ): sample_dir.name
            for sample_dir in sample_dirs
        }
        for future in as_completed(futures):
            name, ok = future.result()
            if not ok:
                failed.append(name)

    if failed:
        print("\n❌ Processing failed for:")
        for n in failed:
            print(f" - {n}")
        raise SystemExit(1)
    else:
        print("\n✨ All samples processed successfully!")


if __name__ == "__main__":
    main()

