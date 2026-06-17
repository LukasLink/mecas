#!/usr/bin/env bash
#SBATCH -J speed                 # chmod +x ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/run_pipeline_sbatch.sh
#SBATCH -A lsteinme                       # sbatch ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/run_pipeline_sbatch.sh
#SBATCH --mem=16g                         #  --dependency=afterok:46151850
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH -t 51:00:00
#SBATCH --qos normal
#SBATCH --array=0
#SBATCH -o /g/steinmetz/link/logs/log_%x_%A_%a.out
#SBATCH -e /g/steinmetz/link/logs/log_%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lukas.link@embl.de

module load R/4.5.2-gfbf-2025b

export R_LIBS_USER=/g/steinmetz/link/R-libs/x86_64-pc-linux-gnu/4.5.2
export PATH="/home/link/.conda/envs/py313/bin:$PATH"

Rscript --vanilla ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/R/cli_count.R ~/Amplicon_barcode_analysis/Lukas_Pipeline/binned_PCR_amplicon_UMI_analysis/config.yaml