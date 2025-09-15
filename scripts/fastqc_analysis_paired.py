#!/usr/bin/env python3

import os
import sys
import subprocess
import logging
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm  # Progress bar

# Toggle development mode
test = False
if test:
    print('Running fastqc_analysis_paired.py in development mode!!')
    input_sample_dir = '/home/jpereira/metaGEM_base/toy_set/dataset/'
    param_threads = 6
    output_fastqc_results_dir = '/home/jpereira/shotgun_snake/Results_test/Data/Fastqc/Fastqc_Analysis'
    output_multiqc_results_dir = '/home/jpereira/shotgun_snake/Results_test/Data/Fastqc/Multiqc_Analysis'
    output_manifest_tsv = '/home/jpereira/shotgun_snake/Results_test/Data/Fastqc/manifest.tsv'
    output_log_file = '/home/jpereira/shotgun_snake/Results_test/Data/Fastqc/fastqc_multiqc_analysis.log'
else:
    # Ensure exactly 6 arguments are provided
    if len(sys.argv) != 7:
        print("Usage: python fastqc_analysis_paired.py "
              "<in_dir> <threads> <fastqc_out> <multiqc_out> <manifest> <log>")
        sys.exit(1)

    input_sample_dir    = sys.argv[1]
    param_threads       = int(sys.argv[2])
    output_fastqc_results_dir  = sys.argv[3]
    output_multiqc_results_dir = sys.argv[4]
    output_manifest_tsv = sys.argv[5]
    output_log_file     = sys.argv[6]

# Echo parameters
print(f'input_sample_dir:           {input_sample_dir}')
print(f'output_fastqc_results_dir:  {output_fastqc_results_dir}')
print(f'output_multiqc_results_dir: {output_multiqc_results_dir}')
print(f'log_file:                   {output_log_file}')

# Ensure output directories exist
os.makedirs(output_fastqc_results_dir, exist_ok=True)
os.makedirs(output_multiqc_results_dir, exist_ok=True)

# Setup logging
if os.path.exists(output_log_file):
    os.remove(output_log_file)
logging.basicConfig(
    filename=output_log_file,
    filemode='w',
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# Verify input directory
if not os.path.isdir(input_sample_dir):
    logging.error(f"The directory {input_sample_dir} does not exist.")
    raise FileNotFoundError(f"The directory {input_sample_dir} does not exist.")

# Gather all FASTQ file paths
fastq_files = []
for sample in os.listdir(input_sample_dir):
    sample_path = os.path.join(input_sample_dir, sample)
    if os.path.isdir(sample_path):
        for fname in os.listdir(sample_path):
            fpath = os.path.join(sample_path, fname)
            if os.path.isfile(fpath) and fpath.endswith((".fastq", ".fastq.gz", ".fq.gz", ".fq")):
                fastq_files.append(fpath)
    else:
        logging.warning(f"Skipping {sample}: Not a directory.")

# Function to run FastQC
def run_fastqc(file_path):
    try:
        cmd = ["fastqc", file_path, "-o", output_fastqc_results_dir]
        logging.info(f"Running: {' '.join(cmd)}")
        with open(output_log_file, "a") as log:
            subprocess.run(cmd, stdout=log, stderr=log, check=True)
        return os.path.basename(file_path), os.path.abspath(file_path)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running FastQC on {file_path}: {e}")
    except Exception as e:
        logging.error(f"Unexpected error processing {file_path}: {e}")
    return None

# Run FastQC in parallel with a progress bar
sample_id_list = []
absolute_file_path_list = []
with ProcessPoolExecutor(max_workers=param_threads) as executor:
    futures = {executor.submit(run_fastqc, fp): fp for fp in fastq_files}
    for future in tqdm(as_completed(futures),
                       total=len(futures),
                       desc="Running FastQC"):
        res = future.result()
        if res:
            sample_id_list.append(res[0])
            absolute_file_path_list.append(res[1])

# Write manifest file
if sample_id_list:
    manifest_df = pd.DataFrame({
        'sample-id': sample_id_list,
        'absolute-filepath': absolute_file_path_list
    })
    manifest_df.to_csv(output_manifest_tsv, sep='\t', index=False)
    logging.info(f"Manifest file created at {output_manifest_tsv}.")
else:
    logging.warning("No FASTQ files were processed; manifest will not be created.")

logging.info("FastQC analysis completed successfully.")

# Aggregate with MultiQC
print('Executing MultiQC Analysis')
multiqc_command = [
    "multiqc", output_fastqc_results_dir,
    "-o", output_multiqc_results_dir,
    "--filename", "index.html",
    "-f"   # force overwrite if index.html already exists
]
print(os.listdir(output_multiqc_results_dir))
print(f"Running MultiQC command: {' '.join(multiqc_command)}")
logging.info(f"Running MultiQC: {' '.join(multiqc_command)}")

try:
    with open(output_log_file, "a") as log:
        subprocess.run(multiqc_command, stdout=log, stderr=log, check=True)
    logging.info("MultiQC analysis completed successfully.")
except subprocess.CalledProcessError as e:
    logging.error(f"Error while running MultiQC: {e}")
except Exception as e:
    logging.error(f"Unexpected error during MultiQC execution: {e}")
