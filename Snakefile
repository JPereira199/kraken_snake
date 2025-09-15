import os, re
from snakemake.io import glob_wildcards, expand

# Paths
PATH   = os.getcwd()
S_PATH = os.path.join(PATH, "scripts")

RESULTS_P = config["results.Dir"]
os.makedirs(RESULTS_P, exist_ok=True)

D_PATH = os.path.join(RESULTS_P, "Data")
V_PATH = os.path.join(RESULTS_P, "Visuals")
R_PATH = os.path.join(RESULTS_P, "Report")
W_PATH = os.path.join(R_PATH, "Results")

print(f"D_PATH: {D_PATH}")
print(f"V_PATH: {V_PATH}")

for path_dir in [RESULTS_P, D_PATH, V_PATH]:
    os.makedirs(path_dir, exist_ok=True)

# Input files
FASTQ_D       = config["samples.Dir"]
USER_METADATA = config["metadata.File"]

# Rule subdirectories
FASTQC_D    = "FastQC"
METADATA_D  = "Metadata"
FASTP_D     = "Fastp"
KRAKEN2_D   = "kraken2"
BRACKEN_D   = "bracken"
M_BRACKEN_D = "merge_bracken"
DESEQ2_D    = "deseq2"

# Parallelization levels
LEVELS = list(config["bracken.levels"])

print(f"LEVELS: {LEVELS}")

# Main output files of each rule
rule_files_path = {
    "fastqc_analysis.rule"     : (D_PATH, FASTQC_D, "manifest.tsv"),
    "metadata_validation.rule" : (D_PATH, METADATA_D, "metadata.tsv"),
    "fastp.rule"               : (D_PATH, FASTP_D, ".sentinel.fastp"),
    "kraken2.rule"             : (D_PATH, KRAKEN2_D, ".sentinel.kraken2"),
    "bracken.rule"             : [(D_PATH, BRACKEN_D, f"{level}", ".sentinel.bracken") for level in LEVELS],
    "merge_bracken.rule"       : [(D_PATH, M_BRACKEN_D, f"{level}", ".sentinel.merge_bracken") for level in LEVELS],
    "deseq2.rule"              : [(D_PATH, DESEQ2_D, f"{level}", ".sentinel.deseq2") for level in LEVELS]
}

print(rule_files_path)

final_files = []
for rule_name, pattern in rule_files_path.items():
    if config.get(rule_name, False):
        if isinstance(pattern, list):  # for bracken/merge_bracken (one per level)
            for path_tuple in pattern:
                final_files.append(os.path.join(*path_tuple))
        else:  # for all other rules
            final_files.append(os.path.join(*pattern))

#  Print all final files
print("\nPrinting final files:")
for file in final_files:
    print(file)

# Define the 'all' rule with final_files as input
rule all:
    input:
        final_files

rule fastqc_analysis_paired:
    input:
        fastqc_analysis = os.path.join(S_PATH, 'fastqc_analysis_paired.py'),
        FASTQ_D = FASTQ_D
    threads:
        config['fastqc_analysis_paired.threads']
    output:
        fastqc_dir  = directory(os.path.join(D_PATH,  FASTQC_D, 'FastQC_results')),
        multiqc_dir  = directory(os.path.join(R_PATH,  FASTQC_D, 'MultiQC_results')),
        #multiqc_report = os.path.join(R_PATH,  FASTQC_D, 'MultiQC_results' ,'index.html'), ## In some symstems multiqc does not work (Illegal Construction core dumped)
        Manifest      = os.path.join(D_PATH,   FASTQC_D, "manifest.tsv"),
        log_file      = os.path.join(D_PATH,   FASTQC_D, "fastqc_multiqc_analysis.log")
    conda:
        "envs/qc_env.yaml"
    shell:
        """
        python {input.fastqc_analysis} \
            {input.FASTQ_D} \
            {threads} \
            {output.fastqc_dir} \
            {output.multiqc_dir} \
            {output.Manifest} \
            {output.log_file}
        """
print('metadata_validation rule inputs:')
print(USER_METADATA)
print(FASTQ_D)
print(os.path.join(D_PATH,   FASTQC_D, "manifest.tsv"))
rule metadata_validation:
    input:
        metadata_validation_py = os.path.join(S_PATH, 'metadata_validation_paired.py'),
        metadata = USER_METADATA,
        samples_dir = FASTQ_D,  # Renamed for clarity
        manifest = os.path.join(D_PATH,  FASTQC_D, "manifest.tsv")
    params:
        file_separator = config['metadata_validation.file_separator'],
        sample_id_column = config['metadata_validation.sample_id_column'],
        fill_nan_values = config['metadata_validation.fill_nan_values'],
        paired = config['metadata_validation.paired']
    output:
        validated_metadata = os.path.join(D_PATH,  METADATA_D, 'metadata.tsv'),
        report_metadata = os.path.join(R_PATH,  METADATA_D, 'index.html'),
        output_logs = os.path.join(D_PATH,  METADATA_D, 'log.metadata_validation.txt')
    conda:
        "envs/py_vis_env.yaml"
    shell:
        """
        python {input.metadata_validation_py} \
            --input_metadata_path {input.metadata} \
            --input_samples_directory {input.samples_dir} \
            --input_manifest_path {input.manifest} \
            --param_user_separator "{params.file_separator}" \
            --param_sample_identifier "{params.sample_id_column}" \
            --param_fill_nan_values "{params.fill_nan_values}" \
            --param_paired {params.paired} \
            --output_validated_metadata_tsv {output.validated_metadata} \
            --output_report_metadata_html {output.report_metadata} \
            --output_logs {output.output_logs} 
        """

rule fastp:
    input:
        fastp_py = os.path.join(S_PATH, 'fastp_subsample_parallel_jobs.py'),
        raw_dir   = FASTQ_D
    params:
        jobs = config['fastp.jobs'],
        subsample = config['fastp.subsample']
    threads:
        config['fastp.threads_fastp']
    output:
        fastp_dir      = directory(os.path.join(D_PATH, "Fastp")),
        sentinel_fastp = os.path.join(D_PATH, "Fastp", ".sentinel.fastp")
    log:
        "logs/fastp.log"
    conda:
        "envs/fastp_env.yaml"
    shell:
        """
        python3 {input.fastp_py} \
            --input-raw-dir {input.raw_dir} \
            --params-subsample {params.subsample} \
            --output-dir {output.fastp_dir} \
            --threads {threads} \
            > {log} 2>&1

        touch {output.sentinel_fastp} 
        """

rule kraken2:
    input:
        kraken2_py = os.path.join(S_PATH, 'kraken2_paired.py'),
        raw_dir    = rules.fastp.output.fastp_dir
    params:
        kraken2_db = config['kraken2.kraken2_db'],
        memory_mapping = config['kraken2.memory_mapping']
    threads:
        config['kraken2.threads']
    output:
        kraken2_dir       = directory(os.path.join(D_PATH, KRAKEN2_D)),
        sentinel_kraken2  = os.path.join(D_PATH, KRAKEN2_D, '.sentinel.kraken2')
    log:
        "logs/kraken2.log"
    conda:
        "envs/kraken2_env.yaml"
    shell:
        """
        python3 {input.kraken2_py} \
            --input-raw-dir {input.raw_dir} \
            --params-kraken2-db {params.kraken2_db} \
            --params-memory-mapping {params.memory_mapping} \
            --params-threads {threads} \
            --output-dir {output.kraken2_dir}

        touch {output.sentinel_kraken2}
        """

rule bracken:
    input:
        bracken_py = os.path.join(S_PATH, "bracken_estimation.py"),
        kraken_dir = rules.kraken2.output.kraken2_dir
    params:
        bracken_db     = config["bracken.bracken_db"],
        read_length    = config["bracken.read_length"],
        kraken_db_kmer = config["bracken.kraken_db_kmer"],
    threads:
        config["bracken.threads"]
    output:
        # separate output directory (and sentinel) per level
        bracken_dir      = directory(os.path.join(D_PATH, BRACKEN_D, "{level}")),
        sentinel_bracken = os.path.join(D_PATH, BRACKEN_D, "{level}", ".sentinel.bracken")
    log:
        "logs/bracken_{level}.log"
    conda:
        "envs/bracken_env.yaml"
    shell:
        """
        python3 {input.bracken_py} \
            --input-kraken-dir {input.kraken_dir} \
            --bracken-db {params.bracken_db} \
            --read-length {params.read_length} \
            --kraken-db-kmer {params.kraken_db_kmer} \
            --level {wildcards.level} \
            --threads {threads} \
            --output-dir {output.bracken_dir} \
        &> {log}
        touch {output.sentinel_bracken}
        """


#rule bracken:
#    input:
#        bracken_py     = os.path.join(S_PATH, 'bracken_estimation.py'),
#        kraken_dir     = rules.kraken2.output.kraken2_dir
#    params:
#        bracken_db     = config['bracken.bracken_db'],
#        read_length    = config['bracken.read_length'],
#        level          = config['bracken.level'],
#        kraken_db_kmer = config['bracken.kraken_db_kmer']
#    threads:
#        config['bracken.threads']
#    output:
#        bracken_dir      = directory(os.path.join(D_PATH, BRACKEN)),
#        sentinel_bracken = os.path.join(D_PATH, BRACKEN, ".sentinel.bracken")
#    log:
#        "logs/bracken.log"
#    conda:
#        "envs/bracken_env.yaml"
#    shell:
#        """
#        python3 {input.bracken_py} \
#         --input-kraken-dir {input.kraken_dir} \
#         --bracken-db {params.bracken_db} \
#         --read-length {params.read_length} \
#         --kraken-db-kmer {params.kraken_db_kmer} \
#         --level {params.level} \
#         --threads {threads} \
#         --output-dir {output.bracken_dir} 
#
#        touch {output.sentinel_bracken}
#        """

metrics_str = " ".join(config['merge_bracken.metrics'])
rule merge_bracken:
    input:
        script      = os.path.join(S_PATH, "merge_bracken_tables.py"),
        bracken_dir = rules.bracken.output.bracken_dir
    params:
        # no need to pass level here
        metrics_str = metrics_str
    output:
        # create a subdirectory per level, then a sentinel inside it
        merge_bracken_dir = directory(os.path.join(D_PATH, M_BRACKEN_D, "{level}")),
        sentinel   = os.path.join(D_PATH, M_BRACKEN_D, "{level}", ".sentinel.merge_bracken")
    log:
        # separate logs can also help debug per level
        "logs/merge_bracken_{level}.log"
    conda:
        "envs/py_vis_env.yaml"
    shell:
        """
        python3 {input.script} \
            --input-dir {input.bracken_dir} \
            --level {wildcards.level} \
            --metrics {params.metrics_str} \
            --output-dir {output.merge_bracken_dir} \
        &> {log}
        touch {output.sentinel}
        """


#rule merge_bracken:
#    input:
#        script        = os.path.join(S_PATH, "merge_bracken_tables.py"),
#        bracken_dir   = rules.bracken.output.bracken_dir
#    params:
#        level         = config["merge_bracken.level"],
#        metrics_str   = metrics_str
#    output:
#        merge_bracken_dir = directory(os.path.join(D_PATH, M_BRACKEN)),
#        sentinel      = os.path.join(D_PATH, M_BRACKEN, ".sentinel.merge_bracken")
#    log:
#        "logs/merge_bracken.log"
#    conda:
#        "envs/py_vis_env.yaml"
#    shell:
#        """
#        python3 {input.script} \
#            --input-dir {input.bracken_dir} \
#            --level {params.level} \
#            --metrics {params.metrics_str} \
#            --output-dir {output.merge_bracken_dir}
#        touch {output.sentinel}
#        """

# TODO: Use internal Paralelization of DESEQ2, given that sometimes fails its excecution
rule deseq2:
    input:
        # Path to your R script
        rscript   = os.path.join(S_PATH, "deseq2_differential_abundance.R"),
        merge_bracken_dir   = rules.merge_bracken.output.merge_bracken_dir,
        metadata  = rules.metadata_validation.output.validated_metadata
    params:
        metric        = config["deseq2.bracken_metric"],
        #level         = config["deseq2.bracken_level"],
        condition     = config["deseq2.condition"],
        batch         = config["deseq2.batch"]
    output:
        deseq2_dir = directory(os.path.join(D_PATH, DESEQ2_D, "{level}")),
        sentinel   = os.path.join(D_PATH, DESEQ2_D, "{level}", ".sentinel.deseq2")
    log:
        "logs/deseq2_{level}.log"
    conda:
        "envs/deseq2_env.yaml"
    shell:
        """
        mkdir -p {output.deseq2_dir}
        awk 'NR != 2' {input.metadata} > {output.deseq2_dir}/deseq2_metadata.tsv

        Rscript {input.rscript} \
            --counts {input.merge_bracken_dir}/bracken_merged_{wildcards.level}_{params.metric}.csv \
            --metadata {output.deseq2_dir}/deseq2_metadata.tsv \
            --condition {params.condition} \
            --batch "{params.batch}" \
            --output {output.deseq2_dir} 
        touch {output.sentinel}
        """

