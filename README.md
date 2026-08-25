# MECAS: Multi-bin Error-Corrected Analysis of Screens

MECAS is a [Snakemake](https://snakemake.readthedocs.io/en/stable/)-based pipeline designed for analysing pooled CRISPR screens of bin sorted samples. 
It processes fastq sequencing files into guide and gene level statistics for overal enrichment trends across bins. 
It uses [bcwithqc](https://github.com/hawkjo/bcwithqc) for error corrected detection of guide sequences and UMIs,
and [MAUDE](https://github.com/de-Boer-Lab/MAUDE?tab=readme-ov-file) for guide and gene level statistical calculations. 


## Installation: 
MECAS works on linux, requires python >= 3.13, and can be installed into your local environment using pip: 
pip install git+https://github.com/LukasLink/mecas.git

### Tips
If you are installing/using MECAS on a cluster you should make sure the environment it is installed to is available to both the login and the compute nodes. 

on first run it will download an apptainer [container](oras://ghcr.io/lukaslink/mecas-container:0.5.0) with all dependencies to:
~/.cache/mecas/apptainer, you can reasign this download location with export MECAS_APPTAINER_PREFIX your/prefered/path

## Usage

### Quick start
Once MECAS is installed you can do a ~10min test run with the provided example files. 