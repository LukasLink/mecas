# Containerization path for the Snakemake pipeline

## 1. Build the image

From the pipeline repository root:

```bash
mkdir -p apptainer containers
cp bcwithqc_maude.def apptainer/bcwithqc_maude.def
apptainer build --fakeroot containers/bcwithqc_maude.sif apptainer/bcwithqc_maude.def
```

If `--fakeroot` is not available on the cluster, build the image somewhere else where Apptainer/Docker builds are allowed, then copy the `.sif` file to `containers/bcwithqc_maude.sif` on the cluster.

## 2. Test the image manually

```bash
apptainer exec containers/bcwithqc_maude.sif Rscript --vanilla -e 'library(tidyverse); library(MAUDE); library(AnnotationDbi); library(org.Hs.eg.db); sessionInfo()'
apptainer exec containers/bcwithqc_maude.sif seqtk 2>&1 | head
apptainer exec containers/bcwithqc_maude.sif STAR --version
apptainer exec containers/bcwithqc_maude.sif bcwithqc --help
```

## 3. Patch the Snakefile

Apply the changes in `Snakefile_container_changes.patch` manually or with `git apply` after copying it into the repo:

```bash
git apply Snakefile_container_changes.patch
```

The key changes are:

- add a global `container: "containers/bcwithqc_maude.sif"`
- remove `module load R/...`
- remove `R_LIBS_USER=...`
- use `Rscript` from inside the container

## 4. Adjust `config.yaml`

For the containerized version, set the bcwithqc executable to the container executable name rather than an absolute host path:

```yaml
bcwithqc:
  bcwithqc_bin: "bcwithqc"
```

Also make sure your QC filtering settings do not ask the shell script to load a `seqtk` module. Inside the container, `seqtk` should be found directly on `$PATH`.

## 5. Run with Snakemake

Dry-run first:

```bash
snakemake -n -p --sdm apptainer
```

Then test only setup:

```bash
snakemake setup -j 1 -p --sdm apptainer \
  --apptainer-args "--bind /g:/g"
```

Then test the whole pipeline:

```bash
snakemake --profile slurm -p --rerun-incomplete --sdm apptainer \
  --apptainer-args "--bind /g:/g"
```

Add further bind mounts if your input/output/config files live outside `/g`, for example `/scratch` or `/data`.
