configfile: "config.yaml"

import csv
import os
import shutil
import math
from pathlib import Path


shell.executable("/usr/bin/bash")

container: "containers/bcwithqc_maude.sif"


#-------------------------------------------------------------------------------
# Merge CLI overrides into nested config
#-------------------------------------------------------------------------------

config.setdefault("paths", {})
config.setdefault("counting", {})

CLI_TO_CONFIG = {
    "output_folder": ("paths", "output_folder"),
    "input_folder": ("paths", "input_folder"),
    "library_path": ("paths", "library_path"),
    "bcwithqc_config_path": ("paths", "bcwithqc_config_path"),
    "fastq_name_table_xlsx": ("paths", "fastq_name_table_xlsx"),
    "data_type": ("counting", "data_type"),
}

for cli_key, (section, key) in CLI_TO_CONFIG.items():
    if cli_key in config:
        config[section][key] = config[cli_key]


#-------------------------------------------------------------------------------
# Validate required user configuration
#-------------------------------------------------------------------------------

required = [
    ("paths", "output_folder"),
    ("paths", "input_folder"),
    ("paths", "library_path"),
    ("paths", "bcwithqc_config_path"),
    ("paths", "fastq_name_table_xlsx"),
    ("counting", "data_type"),
]

missing = []

for section, key in required:
    value = config.get(section, {}).get(key)

    if value is None or str(value).strip() == "":
        missing.append(f"{section}.{key}")

if missing:
    raise ValueError(
        "Missing required configuration values:\n  - "
        + "\n  - ".join(missing)
    )


#-------------------------------------------------------------------------------
# Validate filesystem paths
#-------------------------------------------------------------------------------

def require_file(path, description):
    if not os.path.isfile(path):
        raise ValueError(
            f"{description} does not exist or is not a file:\n{path}"
        )


def require_dir(path, description):
    if not os.path.isdir(path):
        raise ValueError(
            f"{description} does not exist or is not a directory:\n{path}"
        )


require_dir(
    config["paths"]["input_folder"],
    "Input folder"
)

require_file(
    config["paths"]["library_path"],
    "Library file"
)

require_file(
    config["paths"]["bcwithqc_config_path"],
    "bcwithqc config"
)

require_file(
  config["paths"]["fastq_name_table_xlsx"],
    "FASTQ name table"
)


#-------------------------------------------------------------------------------
# Validate data type
#-------------------------------------------------------------------------------

data_type = str(
    config["counting"]["data_type"]
).lower()

if data_type not in {"reads", "umis"}:
    raise ValueError(
        "counting.data_type must be either 'reads' or 'umis'. "
        f"Found '{data_type}'."
    )

config["counting"]["data_type"] = data_type

#-------------------------------------------------------------------------------
# Define paths based on config
#-------------------------------------------------------------------------------
OUTPUT_FOLDER = config["paths"]["output_folder"]
STATE_DIR = os.path.join(OUTPUT_FOLDER, ".pipeline_state")

FASTQ_SYMLINK_FOLDER = os.path.join(OUTPUT_FOLDER, "fastq_symlinks")
BCWITHQC_INPUT_FOLDER = os.path.join(OUTPUT_FOLDER, "bcwithqc_symlinks")
BCWITHQC_OUTPUT_FOLDER = os.path.join(OUTPUT_FOLDER, "bcwithqc_output")
QC_FILTERED_FOLDER = os.path.join(OUTPUT_FOLDER, "QC_filtered")
RDS_OUTPUT_FOLDER = os.path.join(OUTPUT_FOLDER, "rds")
RESULTS_OUTPUT_FOLDER = os.path.join(OUTPUT_FOLDER, "results")
PLOTS_OUTPUT_FOLDER = os.path.join(OUTPUT_FOLDER, "plots")

BCWITHQC_CONFIG = config["paths"]["bcwithqc_config_path"]
BCWITHQC_DATA_TYPE = config["counting"]["data_type"]
#-------------------------------------------------------------------------------
# Manifest helpers
#
# Manifest design:
#   - one row per FASTQ file
#   - pipeline_name identifies the sample/library
#   - read is empty for single-end, or R1/R2 for paired-end
#   - fastq_id uniquely identifies an individual FASTQ file
#-------------------------------------------------------------------------------

BCWITHQC_READ_MATRIX_OUTPUTS = [
    os.path.join(
        BCWITHQC_OUTPUT_FOLDER,
        "{sample}",
        "raw_reads_bc_matrix",
        filename,
    )
    for filename in (
        "barcodes.tsv.gz",
        "features.tsv.gz",
        "matrix.mtx.gz",
    )
]


BCWITHQC_UMI_MATRIX_OUTPUTS = [
    os.path.join(
        BCWITHQC_OUTPUT_FOLDER,
        "{sample}",
        "raw_umis_bc_matrix",
        filename,
    )
    for filename in (
        "barcodes.tsv.gz",
        "features.tsv.gz",
        "matrix.mtx.gz",
    )
]


BCWITHQC_REQUIRED_MATRIX_OUTPUTS = list(
    BCWITHQC_READ_MATRIX_OUTPUTS
)

if BCWITHQC_DATA_TYPE == "umis":
    BCWITHQC_REQUIRED_MATRIX_OUTPUTS.extend(
        BCWITHQC_UMI_MATRIX_OUTPUTS
    )

def read_fastq_manifest():
    manifest_path = str(checkpoints.setup.get().output.manifest)

    with open(manifest_path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    required_columns = {
        "pipeline_name",
        "fastq_id",
        "read",
        "symlink_file_basename",
        "symlink_file",
    }

    if not rows:
        raise ValueError(f"FASTQ manifest is empty: {manifest_path}")

    missing_columns = required_columns.difference(rows[0].keys())
    if missing_columns:
        raise ValueError(
            "FASTQ manifest is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    fastq_ids = [row["fastq_id"] for row in rows]
    if len(fastq_ids) != len(set(fastq_ids)):
        raise ValueError("FASTQ manifest contains duplicate fastq_id values")

    return rows


def manifest_rows_for_sample(sample):
    rows = [
        row for row in read_fastq_manifest()
        if row["pipeline_name"] == sample
    ]

    if not rows:
        raise ValueError(f"No manifest rows found for sample {sample}")

    return rows
def normalize_manifest_read_value(value):
    read_type = (value or "").strip().upper()

    if read_type == "":
        return "SE"

    if read_type not in {"R1", "R2"}:
        raise ValueError(f"Invalid manifest read value: {value!r}. Expected empty, R1, or R2.")

    return read_type


def get_sample_layout(sample):
    rows = manifest_rows_for_sample(sample)
    rows_by_read = {}

    for row in rows:
        read_type = normalize_manifest_read_value(row["read"])

        if read_type in rows_by_read:
            raise ValueError(f"Sample {sample} contains more than one manifest row for read type {read_type}")

        rows_by_read[read_type] = row

    if set(rows_by_read) == {"SE"}:
        return {
            "mode": "SE",
            "rows": rows,
            "SE": rows_by_read["SE"],
        }

    if set(rows_by_read) == {"R1", "R2"}:
        return {
            "mode": "PE",
            "rows": rows,
            "R1": rows_by_read["R1"],
            "R2": rows_by_read["R2"],
        }

    observed = ", ".join(sorted(rows_by_read))

    raise ValueError(
        f"Invalid FASTQ layout for sample {sample}. Found read types: {observed}. "
        "Expected exactly one SE row or exactly one R1 and one R2 row."
    )


def get_se_samples_after_setup(wildcards):
    return sorted(
        sample for sample in get_samples_after_setup(wildcards)
        if get_sample_layout(sample)["mode"] == "SE"
    )


def get_pe_samples_after_setup(wildcards):
    return sorted(
        sample for sample in get_samples_after_setup(wildcards)
        if get_sample_layout(sample)["mode"] == "PE"
    )


def get_manifest_fastq(sample, read_type):
    layout = get_sample_layout(sample)

    if read_type not in layout:
        raise ValueError(f"Sample {sample} does not contain read type {read_type}")

    path = layout[read_type]["symlink_file"]

    if not os.path.exists(path):
        raise ValueError(f"Manifest FASTQ does not exist for sample {sample}, read {read_type}: {path}")

    return path


def get_se_input_fastq(wildcards):
    return get_manifest_fastq(wildcards.sample, "SE")


def get_pe_input_fastq_r1(wildcards):
    return get_manifest_fastq(wildcards.sample, "R1")


def get_pe_input_fastq_r2(wildcards):
    return get_manifest_fastq(wildcards.sample, "R2")


def get_qc_fastq_path(sample, read_type):
    if read_type == "SE":
        return os.path.join(QC_FILTERED_FOLDER, "SE", sample + ".fastq.gz")

    if read_type in {"R1", "R2"}:
        return os.path.join(QC_FILTERED_FOLDER, "PE", sample + "_" + read_type + ".fastq.gz")

    raise ValueError(f"Unsupported read type: {read_type}")


def get_qc_fastqs_for_sample(wildcards):
    layout = get_sample_layout(wildcards.sample)

    return [
        get_qc_fastq_path(wildcards.sample, normalize_manifest_read_value(row["read"]))
        for row in layout["rows"]
    ]


def get_qc_read_count_for_sample(wildcards):
    mode = get_sample_layout(wildcards.sample)["mode"]
    return os.path.join(STATE_DIR, "QC_filter_read_counts", mode, wildcards.sample + ".txt")

def get_samples_after_setup(wildcards):
    samples = sorted({row["pipeline_name"] for row in read_fastq_manifest()})

    if len(samples) == 0:
        raise ValueError(
            "No samples found after setup. Manifest: "
            + str(checkpoints.setup.get().output.manifest)
        )

    return samples


def get_qc_job_done_files(wildcards):
    se_done = expand(
        os.path.join(STATE_DIR, "QC_filter_jobs", "SE", "{sample}.done"),
        sample=get_se_samples_after_setup(wildcards),
    )

    pe_done = expand(
        os.path.join(STATE_DIR, "QC_filter_jobs", "PE", "{sample}.done"),
        sample=get_pe_samples_after_setup(wildcards),
    )

    return se_done + pe_done


def get_prepare_bcwithqc_input_done_files(wildcards):
    return expand(
        os.path.join(STATE_DIR, "bcwithqc_input_done", "{sample}.done"),
        sample=get_samples_after_setup(wildcards),
    )


def get_bcwithqc_done_files(wildcards):
    return expand(
        os.path.join(STATE_DIR, "bcwithqc_done", "{sample}.done"),
        sample=get_samples_after_setup(wildcards),
    )

    
def estimate_qc_filter_runtime(paths):
    input_size_gb = sum(os.path.getsize(str(path)) for path in paths) / (1024 ** 3)

    observed_minutes_per_gb = 0.6
    safety_factor = 1.5
    startup_buffer_minutes = 10

    estimated_runtime = input_size_gb * observed_minutes_per_gb * safety_factor + startup_buffer_minutes

    return math.ceil(estimated_runtime)


def get_qc_filter_se_runtime(wildcards, input):
    return estimate_qc_filter_runtime([input.fastq])


def get_qc_filter_pe_runtime(wildcards, input):
    return estimate_qc_filter_runtime([input.fastq_r1, input.fastq_r2])
  
def get_bcwithqc_preprocess_runtime(wildcards, input):
    with open(str(input.read_count), "r") as handle:
        retained_read_pairs = int(handle.read().strip())

    read_pairs_million = retained_read_pairs / 1_000_000

    estimated_minutes = 60 + 5 * read_pairs_million

    return estimated_minutes


def get_bcwithqc_count_runtime(wildcards, input):
    with open(str(input.read_count), "r") as handle:
        retained_read_pairs = int(handle.read().strip())

    read_pairs_million = retained_read_pairs / 1_000_000

    estimated_minutes = 60 + 25 * read_pairs_million 

    return estimated_minutes
    
#-------------------------------------------------------------------------------
# SNAKEMAKE RULES
#-------------------------------------------------------------------------------
rule all:
    input:
        done_01 = os.path.join(STATE_DIR, "01_setup.done"),
        done_02 = os.path.join(STATE_DIR, "02_infer_QC_filter_params.done"),
        done_03 = os.path.join(STATE_DIR, "03_QC_filter.done"),
        done_04 = os.path.join(STATE_DIR, "04_prepare_bcwithqc_input.done"),
        done_05 = os.path.join(STATE_DIR, "05_bcwithqc.done"),
        done_06 = os.path.join(STATE_DIR, "06_count.done"),
        done_07 = os.path.join(STATE_DIR, "07_MAUDE.done"),
        done_08 = os.path.join(STATE_DIR, "08_plot.done")


checkpoint setup:
    output:
        cfg_rds = os.path.join(STATE_DIR, "resolved_config.rds"),
        cfg_yaml = os.path.join(STATE_DIR, "resolved_config.yaml"),
        manifest = os.path.join(STATE_DIR, "fastq_manifest.tsv"),
        done = os.path.join(STATE_DIR, "01_setup.done")
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 10
    run:
        import yaml

        merged_cfg = os.path.join(STATE_DIR, "merged_config.yaml")
        os.makedirs(STATE_DIR, exist_ok=True)

        with open(merged_cfg, "w") as handle:
            yaml.safe_dump(config, handle, sort_keys=False)

        shell(f"""
            unset R_LIBS
            unset R_LIBS_USER
            Rscript --vanilla R/cli_A_setup.R {merged_cfg}
        """)


rule infer_QC_filter_params:
    input:
        cfg_rds = os.path.join(STATE_DIR, "resolved_config.rds"),
        setup_done = os.path.join(STATE_DIR, "01_setup.done")
    output:
        params = os.path.join(STATE_DIR, "QC_filtering_params.sh"),
        done = os.path.join(STATE_DIR, "02_infer_QC_filter_params.done")
    threads: 1
    resources:
        mem_mb = 1000,
        runtime = 3
    shell:
        r"""
        unset R_LIBS
        unset R_LIBS_USER
        Rscript --vanilla R/cli_B_infer_QC_filter_params.R {input.cfg_rds}
        """


# One QC job per FASTQ file.
#
# Single-end files use QC_MIN_LENGTH_SE.
# Paired-end R1 and R2 files are filtered independently and may use
# different minimum read lengths.
rule QC_filter_SE_worker:
    input:
        fastq = get_se_input_fastq,
        params = os.path.join(STATE_DIR, "QC_filtering_params.sh"),
        setup_done = os.path.join(STATE_DIR, "01_setup.done"),
        infer_QC_filter_params_done = os.path.join(STATE_DIR, "02_infer_QC_filter_params.done")
    output:
        fastq = os.path.join(QC_FILTERED_FOLDER, "SE", "{sample}.fastq.gz"),
        report = os.path.join(STATE_DIR, "QC_filter_reports", "SE", "{sample}.cutadapt.json"),
        read_count = os.path.join(STATE_DIR, "QC_filter_read_counts", "SE", "{sample}.txt"),
        done = os.path.join(STATE_DIR, "QC_filter_jobs", "SE", "{sample}.done")
    threads: 4
    resources:
        mem_mb = 2000,
        runtime = get_qc_filter_se_runtime
    shell:
        r"""
        source "{input.params}"
        
        echo "============================================================"
        echo "JOB RESOURCES"
        echo "Rule:       QC_filter_SE_worker"
        echo "Sample:     {wildcards.sample}"
        echo "Threads:    {threads}"
        echo "Memory:     {resources.mem_mb} MB"
        echo "Runtime:    {resources.runtime} minutes"
        echo "============================================================"

        if [[ "$QC_FILTERING_RUN" == "true" ]]; then
          if [[ -z "$QC_MIN_LENGTH_SE" ]]; then
            echo "ERROR: No QC minimum length was defined for single-end reads" >&2
            exit 1
          fi

          bash shell/QC_filtering_snake.sh \
            --mode SE \
            --input-fastq "{input.fastq}" \
            --output-fastq "{output.fastq}" \
            --min-qual "$QC_MIN_QUAL" \
            --min-length "$QC_MIN_LENGTH_SE" \
            --json-report "{output.report}" \
            --read-count-output "{output.read_count}" \
            --threads "{threads}"
        else
          mkdir -p "$(dirname "{output.fastq}")"
          mkdir -p "$(dirname "{output.report}")"
          mkdir -p "$(dirname "{output.read_count}")"

          count_fastq_records() {{
            local fastq="$1"
            local line_count

            if [[ "$fastq" == *.gz ]]; then
              line_count="$(gzip -cd -- "$fastq" | wc -l)"
            else
              line_count="$(wc -l < "$fastq")"
            fi

            if (( line_count % 4 != 0 )); then
              echo "ERROR: FASTQ line count is not divisible by four: $fastq" >&2
              echo "Observed line count: $line_count" >&2
              exit 1
            fi

            printf '%s\n' "$((line_count / 4))"
          }}

          read_count="$(count_fastq_records "{input.fastq}")"

          printf '%s\n' "$read_count" > "{output.read_count}"

          printf \
            '{{"qc_filtering_run": false, "read_counts": {{"input": %s, "output": %s}}}}\n' \
            "$read_count" \
            "$read_count" \
            > "{output.report}"

          echo "QC filtering disabled."
          echo "Input/retained reads: $read_count"
        fi

        test -e "{output.fastq}"
        test -s "{output.report}"
        test -s "{output.read_count}"

        mkdir -p "$(dirname "{output.done}")"
        touch "{output.done}"
        """
        
rule QC_filter_PE_worker:
    input:
        fastq_r1 = get_pe_input_fastq_r1,
        fastq_r2 = get_pe_input_fastq_r2,
        params = os.path.join(STATE_DIR, "QC_filtering_params.sh"),
        setup_done = os.path.join(STATE_DIR, "01_setup.done"),
        infer_QC_filter_params_done = os.path.join(STATE_DIR, "02_infer_QC_filter_params.done")
    output:
        fastq_r1 = os.path.join(QC_FILTERED_FOLDER, "PE", "{sample}_R1.fastq.gz"),
        fastq_r2 = os.path.join(QC_FILTERED_FOLDER, "PE", "{sample}_R2.fastq.gz"),
        report = os.path.join(STATE_DIR, "QC_filter_reports", "PE", "{sample}.cutadapt.json"),
        read_count = os.path.join(STATE_DIR, "QC_filter_read_counts", "PE", "{sample}.txt"),
        done = os.path.join(STATE_DIR, "QC_filter_jobs", "PE", "{sample}.done")
    threads: 4
    resources:
        mem_mb = 2000,
        runtime = get_qc_filter_pe_runtime
    shell:
        r"""
        source "{input.params}"
        
        echo "============================================================"
        echo "JOB RESOURCES"
        echo "Rule:       QC_filter_PE_worker"
        echo "Sample:     {wildcards.sample}"
        echo "Threads:    {threads}"
        echo "Memory:     {resources.mem_mb} MB"
        echo "Runtime:    {resources.runtime} minutes"
        echo "============================================================"

        if [[ "$QC_FILTERING_RUN" == "true" ]]; then
          if [[ -z "$QC_MIN_LENGTH_R1" || -z "$QC_MIN_LENGTH_R2" ]]; then
            echo "ERROR: QC minimum lengths were not defined for both R1 and R2" >&2
            exit 1
          fi

          bash shell/QC_filtering_snake.sh \
            --mode PE \
            --input-fastq-r1 "{input.fastq_r1}" \
            --input-fastq-r2 "{input.fastq_r2}" \
            --output-fastq-r1 "{output.fastq_r1}" \
            --output-fastq-r2 "{output.fastq_r2}" \
            --min-qual "$QC_MIN_QUAL" \
            --min-length-r1 "$QC_MIN_LENGTH_R1" \
            --min-length-r2 "$QC_MIN_LENGTH_R2" \
            --json-report "{output.report}" \
            --read-count-output "{output.read_count}" \
            --threads "{threads}"
        else
          mkdir -p "$(dirname "{output.fastq_r1}")"
          mkdir -p "$(dirname "{output.report}")"
          mkdir -p "$(dirname "{output.read_count}")"

          if [[ "{input.fastq_r1}" == *.gz ]]; then
            ln -sf "$(realpath "{input.fastq_r1}")" "{output.fastq_r1}"
          else
            gzip -c "{input.fastq_r1}" > "{output.fastq_r1}"
          fi

          if [[ "{input.fastq_r2}" == *.gz ]]; then
            ln -sf "$(realpath "{input.fastq_r2}")" "{output.fastq_r2}"
          else
            gzip -c "{input.fastq_r2}" > "{output.fastq_r2}"
          fi

          count_fastq_records() {{
            local fastq="$1"
            local line_count

            if [[ "$fastq" == *.gz ]]; then
              line_count="$(gzip -cd -- "$fastq" | wc -l)"
            else
              line_count="$(wc -l < "$fastq")"
            fi

            if (( line_count % 4 != 0 )); then
              echo "ERROR: FASTQ line count is not divisible by four: $fastq" >&2
              echo "Observed line count: $line_count" >&2
              exit 1
            fi

            printf '%s\n' "$((line_count / 4))"
          }}

          r1_count="$(count_fastq_records "{input.fastq_r1}")"
          r2_count="$(count_fastq_records "{input.fastq_r2}")"

          if [[ "$r1_count" -ne "$r2_count" ]]; then
            echo "ERROR: Paired FASTQ files contain different numbers of reads" >&2
            echo "R1: $r1_count" >&2
            echo "R2: $r2_count" >&2
            exit 1
          fi

          printf '%s\n' "$r1_count" > "{output.read_count}"

          printf \
            '{{"qc_filtering_run": false, "read_counts": {{"input": %s, "output": %s}}}}\n' \
            "$r1_count" \
            "$r1_count" \
            > "{output.report}"

          echo "QC filtering disabled."
          echo "Input/retained read pairs: $r1_count"
        fi

        test -e "{output.fastq_r1}"
        test -e "{output.fastq_r2}"
        test -s "{output.report}"
        test -s "{output.read_count}"

        mkdir -p "$(dirname "{output.done}")"
        touch "{output.done}"
        """
        
rule QC_filter:
    input:
        get_qc_job_done_files
    output:
        done = os.path.join(STATE_DIR, "03_QC_filter.done")
    threads: 1
    resources:
        mem_mb = 1000,
        runtime = 1
    shell:
        r"""
        touch "{output.done}"
        """

# Collect all QC-filtered FASTQs belonging to one pipeline sample into one
# directory. A paired sample gets both *_R1.fastq.gz and *_R2.fastq.gz there.
rule prepare_bcwithqc_input_worker:
    input:
        manifest = lambda wc: str(checkpoints.setup.get().output.manifest),
        fastqs = get_qc_fastqs_for_sample,
        QC_filter_done = os.path.join(STATE_DIR, "03_QC_filter.done")
    output:
        done = os.path.join(STATE_DIR, "bcwithqc_input_done", "{sample}.done")
    threads: 1
    resources:
        mem_mb = 1000,
        runtime = 3
    run:
        rows = manifest_rows_for_sample(wildcards.sample)
        source_fastqs = list(input.fastqs)

        if len(rows) != len(source_fastqs):
            raise ValueError(
                f"Manifest/input count mismatch for {wildcards.sample}: "
                f"{len(rows)} manifest rows versus {len(source_fastqs)} QC FASTQs"
            )

        sample_dir = os.path.join(BCWITHQC_INPUT_FOLDER, wildcards.sample)

        # Remove stale files, for example when switching a sample from SE to PE.
        if os.path.lexists(sample_dir):
            if os.path.isdir(sample_dir) and not os.path.islink(sample_dir):
                shutil.rmtree(sample_dir)
            else:
                os.unlink(sample_dir)

        os.makedirs(sample_dir, exist_ok=True)

        for row, source_fastq in zip(rows, source_fastqs):
            target_fastq = os.path.join(
                sample_dir,
                row["symlink_file_basename"],
            )
            os.symlink(os.path.realpath(str(source_fastq)), target_fastq)

        Path(os.path.dirname(str(output.done))).mkdir(parents=True, exist_ok=True)
        Path(str(output.done)).touch()
        
rule prepare_bcwithqc_input:
    input:
        get_prepare_bcwithqc_input_done_files
    output:
        done=os.path.join(STATE_DIR, "04_prepare_bcwithqc_input.done")
    threads: 1
    resources:
        mem_mb=1000,
        runtime=1
    shell:
        r"""
        touch "{output.done}"
        """

rule bcwithqc_preprocess:
    input:
        prepared = os.path.join(STATE_DIR, "bcwithqc_input_done", "{sample}.done"),
        preparation_done=os.path.join(STATE_DIR, "04_prepare_bcwithqc_input.done"),
        fastqs = get_qc_fastqs_for_sample,
        read_count = get_qc_read_count_for_sample
    output:
        done = os.path.join(STATE_DIR, "bcwithqc_preprocess_done", "{sample}.done")
    params:
        input_dir = lambda wc: os.path.join(BCWITHQC_INPUT_FOLDER, wc.sample),
        output_dir = lambda wc: os.path.join(BCWITHQC_OUTPUT_FOLDER, wc.sample)
    threads: 10
    resources:
        mem_mb = 10000,
        runtime = get_bcwithqc_preprocess_runtime,
        constraint = "avx512"
    shell:
        r"""
        echo "============================================================"
        echo "JOB RESOURCES"
        echo "Rule:       bcwithqc_preprocess"
        echo "Sample:     {wildcards.sample}"
        echo "Threads:    {threads}"
        echo "Memory:     {resources.mem_mb} MB"
        echo "Runtime:    {resources.runtime} minutes"
        echo "Constraint: {resources.constraint}"
        echo "============================================================"
        
        unset PYTHONPATH
        unset PYTHONHOME
        export PYTHONNOUSERSITE=1

        rm -rf "{params.output_dir}"
        mkdir -p "{params.output_dir}"

        bcwithqc preprocess \
          "{params.input_dir}" \
          --config="{BCWITHQC_CONFIG}" \
          --output-dir="{params.output_dir}" \
          --threads="{threads}" \
          --output-format-bam \
          -vv

        mapfile -t bam_files < <(find "{params.output_dir}" -maxdepth 1 -type f -name "*.bam" | sort)

        if [[ "${{#bam_files[@]}}" -ne 1 ]]; then
          echo "ERROR: Expected exactly one BAM from bcwithqc preprocess, found ${{#bam_files[@]}}" >&2
          echo "Files in {params.output_dir}:" >&2
          find "{params.output_dir}" -maxdepth 2 -type f >&2
          exit 1
        fi
        
        mkdir -p "$(dirname "{output.done}")"
        touch "{output.done}"
        """
        
rule bcwithqc_count:
    input:
        preprocess_done = os.path.join(STATE_DIR, "bcwithqc_preprocess_done", "{sample}.done"),
        read_count = get_qc_read_count_for_sample
    output:
        done = os.path.join(STATE_DIR, "bcwithqc_done", "{sample}.done"),
        matrices = BCWITHQC_REQUIRED_MATRIX_OUTPUTS
    params:
        output_dir = lambda wc: os.path.join(BCWITHQC_OUTPUT_FOLDER, wc.sample)
    threads: 10
    resources:
        mem_mb = 70000,
        runtime = get_bcwithqc_count_runtime,
        constraint = "avx512"
    shell:
        r"""
        echo "============================================================"
        echo "JOB RESOURCES"
        echo "Rule:       bcwithqc_count"
        echo "Sample:     {wildcards.sample}"
        echo "Threads:    {threads}"
        echo "Memory:     {resources.mem_mb} MB"
        echo "Runtime:    {resources.runtime} minutes"
        echo "Constraint: {resources.constraint}"
        echo "============================================================"
        
        unset PYTHONPATH
        unset PYTHONHOME
        export PYTHONNOUSERSITE=1
        
        # Snakemake may create the parent directories of declared outputs.
        # bcwithqc incorrectly interprets their existence as completed matrices.
        rm -rf \
          "{params.output_dir}/raw_reads_bc_matrix" \
          "{params.output_dir}/raw_umis_bc_matrix"
          
        bcwithqc count \
          "{params.output_dir}" \
          --STAR-output-dir="{params.output_dir}" \
          --config="{BCWITHQC_CONFIG}" \
          --output-dir="{params.output_dir}" \
          --threads="{threads}" \
          -vv

        mkdir -p "$(dirname "{output.done}")"
        touch "{output.done}"
        """     
        
rule bcwithqc:
    input:
        get_bcwithqc_done_files
    output:
        done=os.path.join(STATE_DIR, "05_bcwithqc.done")
    threads: 1
    resources:
        mem_mb=1000,
        runtime=1
    shell:
        r"""
        touch "{output.done}"
        """

rule count:
    input:
        cfg_rds = os.path.join(STATE_DIR, "resolved_config.rds"),
        bcwithqc_done = os.path.join(STATE_DIR, "05_bcwithqc.done")
    output:
        count_df = os.path.join(RDS_OUTPUT_FOLDER, "count_df.rds"),
        mapping_results = os.path.join(RESULTS_OUTPUT_FOLDER,"mapping_results.xlsx"),
        done = os.path.join(STATE_DIR, "06_count.done")
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 30
    shell:
        r"""
        unset R_LIBS
        unset R_LIBS_USER
        Rscript --vanilla R/cli_C_count.R {input.cfg_rds}
        """


rule MAUDE:
    input:
        cfg_rds = os.path.join(STATE_DIR, "resolved_config.rds"),
        count_df = os.path.join(RDS_OUTPUT_FOLDER, "count_df.rds"),
        count_done=os.path.join(STATE_DIR, "06_count.done")
    output:
        MAUDE_guide_stats = os.path.join(RDS_OUTPUT_FOLDER, "MAUDE_guide_stats.rds"),
        MAUDE_count_df = os.path.join(RDS_OUTPUT_FOLDER, "MAUDE_ready_count_df.rds"),
        MAUDE_binStats = os.path.join(RDS_OUTPUT_FOLDER, "MAUDE_binStats.rds"),
        MAUDE_gene_stats = os.path.join(RDS_OUTPUT_FOLDER, "MAUDE_gene_stats.rds"),
        done = os.path.join(STATE_DIR, "07_MAUDE.done")
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 120
    shell:
        r"""
        unset R_LIBS
        unset R_LIBS_USER
        Rscript --vanilla R/cli_D_MAUDE.R {input.cfg_rds}
        """


rule plot:
    input:
        cfg_rds = os.path.join(STATE_DIR, "resolved_config.rds"),
        count_df = os.path.join(RDS_OUTPUT_FOLDER, "count_df.rds"),
        MAUDE_count_df = os.path.join(RDS_OUTPUT_FOLDER, "MAUDE_ready_count_df.rds"),
        MAUDE_guide_stats = os.path.join(RDS_OUTPUT_FOLDER, "MAUDE_guide_stats.rds"),
        MAUDE_gene_stats = os.path.join(RDS_OUTPUT_FOLDER, "MAUDE_gene_stats.rds"),
        MAUDE_done = os.path.join(STATE_DIR, "07_MAUDE.done")
    output:
        violin = os.path.join(PLOTS_OUTPUT_FOLDER, "01_read_or_umi_count_plots", "count_summary.xlsx"),
        bar = os.path.join(PLOTS_OUTPUT_FOLDER, "02_MAUDE_QC_plots", "sum_per_sample.png"),
        done = os.path.join(STATE_DIR, "08_plot.done")
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 120
    shell:
        r"""
        unset R_LIBS
        unset R_LIBS_USER
        Rscript --vanilla R/cli_E_plot.R {input.cfg_rds}
        """
