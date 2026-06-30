configfile: "config.yaml"

import os

#-------------------------------------------------------------------------------
# Define paths based on config
#-------------------------------------------------------------------------------  
shell.executable("/usr/bin/bash")

OUTPUT_FOLDER = config["paths"]["output_folder"]
STATE_DIR = os.path.join(OUTPUT_FOLDER, ".pipeline_state")

FASTQ_SYMLINK_FOLDER = os.path.join(OUTPUT_FOLDER, "fastq_symlinks")
BCWITHQC_INPUT_FOLDER = os.path.join(OUTPUT_FOLDER, "bcwithqc_symlinks")
GENOME_OUTPUT_FOLDER = os.path.join(OUTPUT_FOLDER, "ref")
DEDUP_OUTPUT_FOLDER = os.path.join(OUTPUT_FOLDER, "dedup")
BCWITHQC_OUTPUT_FOLDER = os.path.join(OUTPUT_FOLDER, "bcwithqc_output")
MAPPED_OUTPUT_FOLDER = os.path.join(OUTPUT_FOLDER, "mapped")
QC_FILTERED_FOLDER = os.path.join(OUTPUT_FOLDER, "QC_filtered")
RDS_OUTPUT_FOLDER = os.path.join(OUTPUT_FOLDER, "rds")
RESULTS_OUTPUT_FOLDER = os.path.join(OUTPUT_FOLDER, "results")
PLOTS_OUTPUT_FOLDER = os.path.join(OUTPUT_FOLDER, "plots")
#-------------------------------------------------------------------------------
# Helper functions
#-------------------------------------------------------------------------------  
def find_symlink_fastq(sample):
    candidates = [
        os.path.join(FASTQ_SYMLINK_FOLDER, sample + ".fastq.gz"),
        os.path.join(FASTQ_SYMLINK_FOLDER, sample + ".fastq"),
    ]

    matches = [p for p in candidates if os.path.exists(p)]

    if len(matches) == 0:
        raise ValueError(f"No symlink FASTQ found for sample {sample}")
    if len(matches) > 1:
        raise ValueError(f"Both .fastq and .fastq.gz exist for sample {sample}")

    return matches[0]
  
def get_samples_after_setup(wildcards):
    checkpoints.setup.get()
    
    samples = glob_wildcards(
        os.path.join(FASTQ_SYMLINK_FOLDER, "{sample}.fastq.gz")
    ).sample
    
    samples += glob_wildcards(
        os.path.join(FASTQ_SYMLINK_FOLDER, "{sample}.fastq")
    ).sample
    
    samples = sorted(set(samples))
    
    if len(samples) == 0:
      raise ValueError(
          "No samples found after setup_R. Expected symlinks in: "
          + FASTQ_SYMLINK_FOLDER
          + "\nInput folder: "
          + config["paths"]["input_folder"]
          + "\nFASTQ name table: "
          + str(config["paths"].get("fastq_name_table_xlsx", "<not configured>"))
      )
    return samples  
  
def get_bcwithqc_done_files(wildcards):
    samples = get_samples_after_setup(wildcards)
    return expand(
        os.path.join(STATE_DIR, "bcwithqc_done", "{sample}.done"),
        sample=samples
    )
#-------------------------------------------------------------------------------
# SNAKEMAKE RULES
#-------------------------------------------------------------------------------      
rule all:
    input:
        count_df = os.path.join(RDS_OUTPUT_FOLDER, "count_df.rds")

checkpoint setup:
    input:
        config = "config.yaml"
    output:
        cfg_rds = os.path.join(STATE_DIR, "resolved_config.rds"),
        cfg_yaml = os.path.join(STATE_DIR, "resolved_config.yaml"),
        manifest = os.path.join(STATE_DIR, "fastq_manifest.tsv"),
        setup_done = os.path.join(STATE_DIR, "01_setup.done")
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 10
    shell:
        r"""
        module load R/4.5.2-gfbf-2025b
        # export R_LIBS=/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/jupyterhub-ubuntu-24.04/4.5.2:/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/jupyterhub-ubuntu-22.04/4.4.1
        export R_LIBS_USER=/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/4.5.2
        Rscript --vanilla R/cli_setup.R {input.config}
        """

# I AM DEPRECATING THE ENTIRE ALIGN UMI TOOLS BRANCH        
# rule make_reference_R:
#     input:
#         cfg_rds = os.path.join(STATE_DIR, "resolved_config.rds"),
#         setup_done = os.path.join(STATE_DIR, "01_setup.done")
#     output:
#         fasta = os.path.join(GENOME_OUTPUT_FOLDER, "ref.fa"),
#         gtf = os.path.join(GENOME_OUTPUT_FOLDER, "ref.gtf"),
#         STAR_ref_params = os.path.join(STATE_DIR, "STAR_ref_params.sh")
#     threads: 1
#     resources:
#         mem_mb = 4000,
#         runtime = 20
#     shell:
#         r"""
#         module load R/4.5.2-gfbf-2025b
#         # export R_LIBS=/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/jupyterhub-ubuntu-24.04/4.5.2:/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/jupyterhub-ubuntu-22.04/4.4.1
#         export R_LIBS_USER=/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/4.5.2
#         Rscript --vanilla R/cli_make_reference.R {input.cfg_rds} --snakemake
#         """
#         
# rule make_reference_shell:
#     input:
#         fasta = os.path.join(GENOME_OUTPUT_FOLDER, "ref.fa"),
#         gtf = os.path.join(GENOME_OUTPUT_FOLDER, "ref.gtf"),
#         STAR_ref_params = os.path.join(STATE_DIR, "STAR_ref_params.sh")
#     output:
#         done = os.path.join(STATE_DIR, "02_make_reference.done")
#     threads: 12
#     resources:
#         mem_mb = 128000,
#         runtime = 5760
#     shell:
#         r"""
#         module load STAR/2.7.11b
# 
#         source {input.STAR_ref_params}
# 
#         LIMIT_GENOME_GENERATE_RAM=$(( {resources.mem_mb} * 1000000 * 90 / 100 ))
# 
#         mkdir -p "$STAR_INDEX_OUTPUT_DIR"
# 
#         STAR \
#           --runThreadN {threads} \
#           --runMode genomeGenerate \
#           --genomeDir "$STAR_INDEX_OUTPUT_DIR" \
#           --genomeFastaFiles "$REF_FASTA_PATH" \
#           --sjdbGTFfile "$REF_GTF_PATH" \
#           --limitGenomeGenerateRAM "$LIMIT_GENOME_GENERATE_RAM" \
#           --genomeSAindexNbases "$GENOME_SA_INDEX_N_BASES" \
#           --sjdbOverhang "$SJDB_OVERHANG" \
#           --genomeChrBinNbits "$GENOME_CHR_BIN_N_BITS"
# 
#         touch {output.done}
#         """
        
rule QC_filtering_R:
    input:
        cfg_rds = os.path.join(STATE_DIR, "resolved_config.rds"),
        setup_done = os.path.join(STATE_DIR, "01_setup.done")
    output:
        params = os.path.join(STATE_DIR, "QC_filtering_params.sh")
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 10
    shell:
        r"""
        module load R/4.5.2-gfbf-2025b
        export R_LIBS_USER=/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/4.5.2
        Rscript --vanilla R/cli_QC_filtering.R {input.cfg_rds}
        """
        
rule QC_filter_shell:
    input:
        fastq = lambda wc: find_symlink_fastq(wc.sample),
        params = os.path.join(STATE_DIR, "QC_filtering_params.sh"),
        setup_done = os.path.join(STATE_DIR, "01_setup.done")
    output:
        fastq = os.path.join(QC_FILTERED_FOLDER, "{sample}.fastq.gz")
    shell:
        r"""
        source {input.params}

        if [[ "$QC_FILTERING_RUN" == "true" ]]; then
          bash shell/QC_filtering_snake.sh \
            --input-fastq {input.fastq} \
            --output-fastq {output.fastq} \
            --min-qual "$QC_MIN_QUAL" \
            --qual-offset "$QC_QUAL_OFFSET" \
            --min-length "$QC_MIN_LENGTH" \
            --use-modules "$USE_MODULES" \
            --seqtk-module "$SEQTK_MODULE"
        else
          mkdir -p "$(dirname "{output.fastq}")"

          if [[ "{input.fastq}" == *.gz ]]; then
            ln -sf "{input.fastq}" "{output.fastq}"
          else
            gzip -c "{input.fastq}" > "{output.fastq}"
          fi
        fi
        """

rule prepare_bcwithqc_input:
    input:
        fastq = os.path.join(QC_FILTERED_FOLDER, "{sample}.fastq.gz")
    output:
        fastq = os.path.join(BCWITHQC_INPUT_FOLDER, "{sample}", "{sample}.fastq.gz")
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 10
    shell:
        r"""
        mkdir -p "$(dirname {output.fastq})"
        ln -sf "$(realpath "{input.fastq}")" "{output.fastq}"
        """

BCWITHQC_BIN = config["bcwithqc"]["bcwithqc_bin"]
BCWITHQC_CONFIG = config["bcwithqc"]["bcwithqc_config_path"]
        
rule bcwithqc_one_sample:
    input:
        fastq = os.path.join(BCWITHQC_INPUT_FOLDER, "{sample}", "{sample}.fastq.gz")
    output:
        done = os.path.join(STATE_DIR, "bcwithqc_done", "{sample}.done")
    params:
        input_dir = lambda wc: os.path.join(BCWITHQC_INPUT_FOLDER, wc.sample),
        output_dir = lambda wc: os.path.join(BCWITHQC_OUTPUT_FOLDER, wc.sample)
    threads: 10
    resources:
        mem_mb = 40000,
        runtime = 2400
    shell:
        r"""
        unset PYTHONPATH
        unset PYTHONHOME
        export PYTHONNOUSERSITE=1

        rm -rf "{params.output_dir}"
        mkdir -p "{params.output_dir}"

        "{BCWITHQC_BIN}" preprocess \
          "{params.input_dir}" \
          --config="{BCWITHQC_CONFIG}" \
          --output-dir="{params.output_dir}" \
          --threads="{threads}" \
          --output-format-bam \
          -vv

        "{BCWITHQC_BIN}" count \
          "{params.output_dir}" \
          --STAR-output-dir="{params.output_dir}" \
          --config="{BCWITHQC_CONFIG}" \
          --output-dir="{params.output_dir}" \
          --threads="{threads}" \
          -vv
          
        mkdir -p "$(dirname "{output.done}")"
        touch "{output.done}"
        """
        
rule bcwithqc_done:
    input:
        get_bcwithqc_done_files
    output:
        done = os.path.join(STATE_DIR, "05_bcwithqc.done")
    threads: 1
    resources:
        mem_mb = 500,
        runtime = 1
    shell:
        r"""
        touch "{output.done}"
        """
        
rule MAUDE:
    input:
        cfg_rds = os.path.join(STATE_DIR, "resolved_config.rds"),
        setup_done = os.path.join(STATE_DIR, "01_setup.done"),
        bcwithqc_done = os.path.join(STATE_DIR, "05_bcwithqc.done")
    output:
        count_df = os.path.join(RDS_OUTPUT_FOLDER, "count_df.rds"),
        mapping_results = os.path.join(RESULTS_OUTPUT_FOLDER,"mapping_results.xlsx"),
        MAUDE_guide_stats = os.path.join(RDS_OUTPUT_FOLDER, "MAUDE_guide_stats.rds"),
        MAUDE_gene_stats = os.path.join(RDS_OUTPUT_FOLDER, "MAUDE_gene_stats.rds")
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 120
    shell:
        r"""
        module load R/4.5.2-gfbf-2025b
        export R_LIBS_USER=/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/4.5.2
        Rscript --vanilla R/cli_MAUDE.R {input.cfg_rds}
        """
        
rule plot:
    input:
        cfg_rds = os.path.join(STATE_DIR, "resolved_config.rds"),
        setup_done = os.path.join(STATE_DIR, "01_setup.done"),
        count_df = os.path.join(RDS_OUTPUT_FOLDER, "count_df.rds"),
        MAUDE_guide_stats = os.path.join(RDS_OUTPUT_FOLDER, "MAUDE_guide_stats.rds"),
        MAUDE_gene_stats = os.path.join(RDS_OUTPUT_FOLDER, "MAUDE_gene_stats.rds")
    output:
        violin = os.path.join(PLOTS_OUTPUT_FOLDER, "01_read_or_umi_count_plots", "count_summary.xlsx"),
        bar = os.path.join(PLOTS_OUTPUT_FOLDER, "02_MAUDE_QC_plots", "sum_per_sample.png"),
        done = os.path.join(STATE_DIR, "05_plot.done")
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 120
    shell:
        r"""
        module load R/4.5.2-gfbf-2025b
        export R_LIBS_USER=/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/4.5.2
        Rscript --vanilla R/cli_plot.R {input.cfg_rds}
        touch "{output.done}"
        """        
#         
# rule prepare_count_inputs_R:
#   input:
#       cfg_rds = os.path.join(STATE_DIR, "resolved_config.rds"),
#       qc_done = os.path.join(STATE_DIR, "03_QC_filtering.done")
#   output:
#       done = os.path.join(STATE_DIR, "04_prepare_count_inputs.done")
#   shell:
#       r"""
#       module load R/4.5.2-gfbf-2025b
#       export R_LIBS_USER=/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/4.5.2
#       Rscript --vanilla R/cli_prepare_count_inputs.R {input.cfg_rds}
#       touch {output.done}
#       """
