#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path
import argparse
from itertools import chain
from concurrent.futures import ProcessPoolExecutor, as_completed
import tempfile
import shutil

def run_command(cmd):
    """Execute a command and check for success."""
    print(f" ▶︎ {' '.join(cmd)}")
    result = subprocess.run(cmd, check=True)
    return result.returncode

def process_sample(sample_dir: Path, output_dir: Path, subsample: int, threads: int):
    """
    Trim with fastp and subsample paired FASTQ files in a temp folder,
    then move results (and fastp logs) into output_dir/sample_name.
    """
    sample_name = sample_dir.name
    print(f"\n🔬 Sample: {sample_name}")

    # find R1/R2
    r1_list = list(chain(sample_dir.glob("*1.fq.gz"), sample_dir.glob("*1.fastq.gz")))
    r2_list = list(chain(sample_dir.glob("*2.fq.gz"), sample_dir.glob("*2.fastq.gz")))
    if len(r1_list) != 1 or len(r2_list) != 1:
        print(f"⚠️  Could not find exactly one R1/R2 in {sample_dir}")
        return (sample_name, False)
    r1, r2 = r1_list[0], r2_list[0]

    # prepare final output folder
    sample_out = output_dir / sample_name
    sample_out.mkdir(parents=True, exist_ok=True)

    # we'll name the final outputs:
    final_r1 = sample_out / f"{sample_name}_1.subsampled.fq.gz"
    final_r2 = sample_out / f"{sample_name}_2.subsampled.fq.gz"
    fastp_html = sample_out / f"{sample_name}.fastp.html"
    fastp_json = sample_out / f"{sample_name}.fastp.json"

    try:
        # work in a temp dir
        with tempfile.TemporaryDirectory() as tmpstr:
            tmpdir = Path(tmpstr)

            # copy raw reads in
            tmp_r1 = tmpdir / r1.name
            tmp_r2 = tmpdir / r2.name
            shutil.copy2(r1, tmp_r1)
            shutil.copy2(r2, tmp_r2)

            # 1) trim with fastp
            trimmed_r1 = tmpdir / f"{sample_name}_1.trimmed.fq.gz"
            trimmed_r2 = tmpdir / f"{sample_name}_2.trimmed.fq.gz"
            fastp_cmd = [
                "fastp",
                "-i", str(tmp_r1),
                "-I", str(tmp_r2),
                "-o", str(trimmed_r1),
                "-O", str(trimmed_r2),
                "-h", str(fastp_html),
                "-j", str(fastp_json),
                "-w", str(threads)
            ]
            run_command(fastp_cmd)

            # 2) subsample with reformat.sh
            subs_r1 = tmpdir / f"{sample_name}_1.subsampled.fq.gz"
            subs_r2 = tmpdir / f"{sample_name}_2.subsampled.fq.gz"
            reformat_cmd = [
                "reformat.sh",
                f"in1={trimmed_r1}",
                f"in2={trimmed_r2}",
                f"out1={subs_r1}",
                f"out2={subs_r2}",
                f"reads={subsample}",
                f"threads={threads}"
            ]
            run_command(reformat_cmd)

            # move only the final subsampled files into place
            shutil.move(str(subs_r1), str(final_r1))
            shutil.move(str(subs_r2), str(final_r2))

        print(f"✅ Done: {sample_name}")
        return (sample_name, True)

    except subprocess.CalledProcessError as e:
        print(f"❌ Error in sample {sample_name}: {e}")
        return (sample_name, False)

def main():
    parser = argparse.ArgumentParser(
        description="Trim then subsample paired-end FASTQs (fastp + BBMap reformat.sh)"
    )
    parser.add_argument(
        "--input-raw-dir", required=True, type=Path,
        help="Directory containing one subfolder per sample"
    )
    parser.add_argument(
        "--params-subsample", required=True, type=int,
        help="Number of read-pairs to keep per sample"
    )
    parser.add_argument(
        "--output-dir", required=True, type=Path,
        help="Where to write per-sample output folders"
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Threads per sample (clamped to 1–16; default=4)"
    )
    args = parser.parse_args()

    # Count number of sample folders (directories) in the input directory
    input_dir = Path(args.input_raw_dir)
    sample_dirs = [d for d in input_dir.iterdir() if d.is_dir()]
    num_samples = len(sample_dirs)
    threads = args.threads 

    # Calculate jobs based on available CPUs and per-job threads, but do not exceed number of samples
    threads_per_job = max(min(threads // num_samples, 16),1)
    jobs = min(threads // threads_per_job, num_samples)

    print(f"⚙️  Running up to {jobs} samples in parallel × {threads_per_job} threads each (CPUs: {threads})")

    # validate dirs
    if not args.input_raw_dir.exists():
        print(f"⚠️  Input not found: {args.input_raw_dir}")
        raise SystemExit(1)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sample_dirs = sorted(d for d in args.input_raw_dir.iterdir() if d.is_dir())
    if not sample_dirs:
        print(f"⚠️  No sample subdirectories in: {args.input_raw_dir}")
        raise SystemExit(1)

    # execute
    failed = []
    with ProcessPoolExecutor(max_workers=jobs) as executor:
        futures = {
            executor.submit(
                process_sample,
                sample_dir,
                args.output_dir,
                args.params_subsample,
                threads_per_job
            ): sample_dir.name
            for sample_dir in sample_dirs
        }
        for future in as_completed(futures):
            name, ok = future.result()
            if not ok:
                failed.append(name)

    if failed:
        print("\n❌ Subsampling failed for:")
        for n in failed:
            print("  -", n)
        raise SystemExit(1)
    else:
        print("\n✨ All samples completed successfully!")

if __name__ == "__main__":
    main()
