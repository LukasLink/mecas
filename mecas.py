#!/usr/bin/env python3

import argparse
import importlib.resources
import os
from os import path
import shutil
import sys
from pathlib import Path


PIPELINE_ITEMS = [
    "Snakefile",
    "R",
    "shell",
    "containers",
    "sample_profiles",
    "data",
]


def ensure_symlink(src: Path, dst: Path):
    """
    Ensure dst is a symlink pointing to src.

    Existing incorrect links or copied files/directories are replaced.
    Correct links are left untouched.
    """

    if dst.is_symlink():
        if dst.resolve() == src.resolve():
            return
        dst.unlink()

    elif dst.exists():
        print(f"Replacing existing {dst}")

        if dst.is_dir():
            shutil.rmtree(dst)
        else:
            dst.unlink()

    os.symlink(
        src,
        dst,
        target_is_directory=src.is_dir(),
    )


def prepare_pipeline(source_root: Path, pipeline_dir: Path):
    """
    Ensure the pipeline working directory exists.
    """

    pipeline_dir.mkdir(parents=True, exist_ok=True)

    #
    # Give every run its own editable config.yaml.
    # Do not overwrite an existing config.
    #
    config_src = source_root / "config.yaml"
    config_dst = pipeline_dir / "config.yaml"

    if not config_dst.exists():
        shutil.copy2(config_src, config_dst)

    #
    # Ensure all pipeline components are present.
    #
    for item in PIPELINE_ITEMS:

        src = source_root / item
        dest = pipeline_dir / item
        
        if not src.exists():
            sys.exit(
                f"\nERROR\n"
                f"Missing pipeline component:\n"
                f"{src}"
            )

        ensure_symlink(src, dest)
        
def resolve_user_paths(args):
    """
    Resolve all user-supplied paths to absolute canonical paths.
    """

    args.input_dir = str(Path(args.input_dir).resolve())
    args.output_dir = str(Path(args.output_dir).resolve())
    args.library_file = str(Path(args.library_file).resolve())
    args.bcwithqc_config_file = str(Path(args.bcwithqc_config_file).resolve())
    args.experiment_info_file = str(Path(args.experiment_info_file).resolve())


def get_mount_root(path: str) -> str:
    """
    Return the top-level filesystem component to bind.

    Examples
    --------
    /g/steinmetz/data/file.fastq.gz  -> /g
    /scratch/link/run1              -> /scratch
    /home/user/data                 -> /home
    """

    parts = Path(path).resolve().parts

    if len(parts) < 2:
        return "/"

    return "/" + parts[1]


def build_apptainer_bind_string(args, pipeline_root) -> str:
    """
    Build the Apptainer bind string from all user-supplied paths.
    """

    paths = [
        args.input_dir,
        args.output_dir,
        args.library_file,
        args.bcwithqc_config_file,
        args.experiment_info_file,
        str(pipeline_root),
    ]

    mounts = sorted(
        {
            get_mount_root(path)
            for path in paths
        }
    )

    return ",".join(
        f"{mount}:{mount}"
        for mount in mounts
    )

def main():

    parser = argparse.ArgumentParser(
        prog="mecas",
        description="Run the MECAS amplicon barcode analysis pipeline.",
    )

    parser.add_argument(
        "target",
        nargs="?",
        default="all",
        help="Snakemake target (default: all)",
    )

    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing FASTQ files.",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for the analysis.",
    )

    parser.add_argument(
        "--library-file",
        required=True,
        help="Library annotation file.",
    )

    parser.add_argument(
        "--bcwithqc-config-file",
        required=True,
        help="bcwithqc JSON configuration.",
    )

    parser.add_argument(
        "--experiment-info-file",
        required=True,
        help="Experiment information spreadsheet.",
    )

    parser.add_argument(
        "--reads-or-umis",
        choices=("reads", "umis"),
        default="reads",
        help="Use read counts or UMI counts (default: reads).",
    )
    # Forward all unknown options directly to Snakemake.
    args, snakemake_args = parser.parse_known_args()
    resolve_user_paths(args)
    
    # Locate installed pipeline.
    #
    # pipeline_root = Path(importlib.resources.files("mecas"))
    pipeline_root = Path(__file__).resolve().parent

    output_dir = Path(args.output_dir)
    pipeline_dir = output_dir / "pipeline"

    print(f"Preparing pipeline directory:\n{pipeline_dir}\n")

    prepare_pipeline(
        pipeline_root,
        pipeline_dir,
    )
    
    # Bind all relevant paths to apptainer. 
    apptainer_bind_string = build_apptainer_bind_string(args, pipeline_root)
    
    cmd = [
        "snakemake",
        args.target,
    
        "--sdm",
        "apptainer",
    
        "--apptainer-args",
        f"--bind {apptainer_bind_string}",
    
        "--profile",
        "sample_profiles/slurm",
    
        "--config",

        f"input_folder={args.input_dir}",
        f"output_folder={output_dir}",
        f"library_path={args.library_file}",
        f"bcwithqc_config_path={args.bcwithqc_config_file}",
        f"fastq_name_table_xlsx={args.experiment_info_file}",
        f"data_type={args.reads_or_umis}",
    ]

    cmd.extend(snakemake_args)

    print("Launching Snakemake\n")
    print(" ".join(cmd))
    print()

    os.chdir(pipeline_dir)
    
    print(f"Apptainer bind mounts: {apptainer_bind_string}")
    
    # Replace ourselves with Snakemake.
    os.execvp(
        "snakemake",
        cmd,
    )


if __name__ == "__main__":
    main()
