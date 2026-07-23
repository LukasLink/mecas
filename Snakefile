configfile: "config.yaml"

import csv
import os
import shutil
from pathlib import Path


shell.executable("/usr/bin/bash")

container: "containers/bcwithqc_maude.sif"

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

#-------------------------------------------------------------------------------
# Manifest helpers
#
# Manifest design:
#   - one row per FASTQ file
#   - pipeline_name identifies the sample/library
#   - read is empty for single-end, or R1/R2 for paired-end
#   - fastq_id uniquely identifies an individual FASTQ file
#-------------------------------------------------------------------------------
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


def manifest_row_for_fastq_id(fastq_id):
    matches = [
        row for row in read_fastq_manifest()
        if row["fastq_id"] == fastq_id
    ]

    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one manifest row for fastq_id {fastq_id}, "
            f"found {len(matches)}"
        )

    return matches[0]


def manifest_rows_for_sample(sample):
    rows = [
        row for row in read_fastq_manifest()
        if row["pipeline_name"] == sample
    ]

    if not rows:
        raise ValueError(f"No manifest rows found for sample {sample}")

    return rows

def get_read_type_for_fastq_id(wildcards):
    row = manifest_row_for_fastq_id(wildcards.fastq_id)

    read_type = row["read"].strip().upper()

    # Empty read field means single-end
    if read_type == "":
        return "SE"

    if read_type not in {"R1", "R2"}:
        raise ValueError(
            f"Invalid read value for fastq_id {wildcards.fastq_id}: "
            f"{row['read']!r}. Expected an empty value, R1, or R2."
        )

    return read_type
  
def find_symlink_fastq(wildcards):
    row = manifest_row_for_fastq_id(wildcards.fastq_id)
    path = row["symlink_file"]

    if not os.path.exists(path):
        raise ValueError(
            f"Manifest symlink FASTQ does not exist for {wildcards.fastq_id}: {path}"
        )

    return path


def get_samples_after_setup(wildcards):
    samples = sorted({row["pipeline_name"] for row in read_fastq_manifest()})

    if len(samples) == 0:
        raise ValueError(
            "No samples found after setup. Manifest: "
            + str(checkpoints.setup.get().output.manifest)
        )

    return samples


def get_qc_fastqs_for_sample(wildcards):
    return [
        os.path.join(QC_FILTERED_FOLDER, row["fastq_id"] + ".fastq.gz")
        for row in manifest_rows_for_sample(wildcards.sample)
    ]



def get_fastq_ids_after_setup(wildcards):
    return sorted(row["fastq_id"] for row in read_fastq_manifest())


def get_qc_job_done_files(wildcards):
    return expand(
        os.path.join(STATE_DIR, "QC_filter_jobs", "{fastq_id}.done"),
        fastq_id=get_fastq_ids_after_setup(wildcards),
    )


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
    input:
        config = "config.yaml"
    output:
        cfg_rds = os.path.join(STATE_DIR, "resolved_config.rds"),
        cfg_yaml = os.path.join(STATE_DIR, "resolved_config.yaml"),
        manifest = os.path.join(STATE_DIR, "fastq_manifest.tsv"),
        done = os.path.join(STATE_DIR, "01_setup.done")
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 10
    shell:
        r"""
        unset R_LIBS
        unset R_LIBS_USER
        Rscript --vanilla R/cli_setup.R {input.config}
        """


rule infer_QC_filter_params:
    input:
        cfg_rds = os.path.join(STATE_DIR, "resolved_config.rds"),
        setup_done = os.path.join(STATE_DIR, "01_setup.done")
    output:
        params = os.path.join(STATE_DIR, "QC_filtering_params.sh"),
        done = os.path.join(STATE_DIR, "02_infer_QC_filter_params.done")
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 10
    shell:
        r"""
        unset R_LIBS
        unset R_LIBS_USER
        Rscript --vanilla R/cli_infer_QC_filter_params.R {input.cfg_rds}
        """


# One QC job per FASTQ file.
#
# Single-end files use QC_MIN_LENGTH_SE.
# Paired-end R1 and R2 files are filtered independently and may use
# different minimum read lengths.
rule QC_filter_prep:
    input:
        fastq = find_symlink_fastq,
        params = os.path.join(STATE_DIR, "QC_filtering_params.sh"),
        setup_done = os.path.join(STATE_DIR, "01_setup.done"),
        infer_QC_filter_params_done = os.path.join(STATE_DIR, "02_infer_QC_filter_params.done")
    output:
        fastq = os.path.join(QC_FILTERED_FOLDER, "{fastq_id}.fastq.gz"),
        done = os.path.join(STATE_DIR, "QC_filter_jobs", "{fastq_id}.done")
    params:
        read_type = get_read_type_for_fastq_id
    shell:
        r"""
        source "{input.params}"

        if [[ "$QC_FILTERING_RUN" == "true" ]]; then

          case "{params.read_type}" in
            SE)
              QC_MIN_LENGTH_SELECTED="$QC_MIN_LENGTH_SE"
              ;;
            R1)
              QC_MIN_LENGTH_SELECTED="$QC_MIN_LENGTH_R1"
              ;;
            R2)
              QC_MIN_LENGTH_SELECTED="$QC_MIN_LENGTH_R2"
              ;;
            *)
              echo "ERROR: Unknown read type: {params.read_type}" >&2
              exit 1
              ;;
          esac

          if [[ -z "$QC_MIN_LENGTH_SELECTED" ]]; then
            echo \
              "ERROR: No QC minimum length was defined for read type {params.read_type}" \
              >&2
            exit 1
          fi

          bash shell/QC_filtering_snake.sh \
            --input-fastq "{input.fastq}" \
            --output-fastq "{output.fastq}" \
            --min-qual "$QC_MIN_QUAL" \
            --qual-offset "$QC_QUAL_OFFSET" \
            --min-length "$QC_MIN_LENGTH_SELECTED"

        else
          mkdir -p "$(dirname "{output.fastq}")"

          if [[ "{input.fastq}" == *.gz ]]; then
            ln -sf \
              "$(realpath "{input.fastq}")" \
              "{output.fastq}"
          else
            gzip -c "{input.fastq}" > "{output.fastq}"
          fi
        fi
        
        test -e "{output.fastq}"

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
        mem_mb=1000,
        runtime=1
    shell:
        r"""
        touch "{output.done}"
        """

# Collect all QC-filtered FASTQs belonging to one pipeline sample into one
# directory. A paired sample gets both *_R1.fastq.gz and *_R2.fastq.gz there.
rule prepare_bcwithqc_input_prep:
    input:
        manifest = lambda wc: str(checkpoints.setup.get().output.manifest),
        fastqs = get_qc_fastqs_for_sample,
        QC_filter_done = os.path.join(STATE_DIR, "03_QC_filter.done")
    output:
        done = os.path.join(STATE_DIR, "bcwithqc_input_done", "{sample}.done")
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 10
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

rule bcwithqc_prep:
    input:
        prepared = os.path.join(STATE_DIR, "bcwithqc_input_done", "{sample}.done"),
        preparation_done=os.path.join(STATE_DIR, "04_prepare_bcwithqc_input.done")
    output:
        done = os.path.join(STATE_DIR, "bcwithqc_done", "{sample}.done")
    params:
        input_dir = lambda wc: os.path.join(BCWITHQC_INPUT_FOLDER, wc.sample),
        output_dir = lambda wc: os.path.join(BCWITHQC_OUTPUT_FOLDER, wc.sample)
    threads: 10
    resources:
        # mem_mb = 40000,
        # runtime = 2400,
        mem_mb = 10000,
        runtime = 240,
        constraint = "avx512"
    shell:
        r"""
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

        bam_file="${{bam_files[0]}}"

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
        Rscript --vanilla R/cli_count.R {input.cfg_rds}
        """


rule MAUDE:
    input:
        cfg_rds = os.path.join(STATE_DIR, "resolved_config.rds"),
        count_df = os.path.join(RDS_OUTPUT_FOLDER, "count_df.rds"),
        count_done=os.path.join(STATE_DIR, "06_count.done")
    output:
        MAUDE_guide_stats = os.path.join(RDS_OUTPUT_FOLDER, "MAUDE_guide_stats.rds"),
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
        Rscript --vanilla R/cli_MAUDE.R {input.cfg_rds}
        """


rule plot:
    input:
        cfg_rds = os.path.join(STATE_DIR, "resolved_config.rds"),
        count_df = os.path.join(RDS_OUTPUT_FOLDER, "count_df.rds"),
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
        Rscript --vanilla R/cli_plot.R {input.cfg_rds}
        """
