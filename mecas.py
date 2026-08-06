#!/usr/bin/env python3

import argparse
# import importlib.resources
import os
from os import path
import shutil
import sys
from pathlib import Path
import yaml


VERSION = "0.1.0"


PIPELINE_ITEMS = [
    "Snakefile",
    "R",
    "shell",
    "containers",
    "sample_profiles",
    "data",
]


COMMAND_SPECS = {
    "all": {
        "help": "Run the full pipeline.",
        "description": "Run the complete MECAS pipeline from setup through plotting.",
        "snakemake_target": "all",
    },
    "setup": {
        "help": "Resolve config and create the manifest.",
        "description": "Resolve configuration, validate inputs, and create the setup outputs and FASTQ manifest.",
        "snakemake_target": "setup",
    },
    "QC_filter": {
        "help": "Run FASTQ QC filtering.",
        "description": "Quality filter the input FASTQ files and record retained read counts.",
        "snakemake_target": "QC_filter",
    },
    "bcwithqc_preprocess": {
        "help": "Prepare input for bcwithqc.",
        "description": "Run bcwithqc preprocess to Error correct the barcodes/guides specified int the bcwithqc config file.",
        "snakemake_target": "bcwithqc_preprocess",
    },
    "bcwithqc_count": {
        "help": "Run bcwithqc counting.",
        "description": "Run bcwithqc count to count reads, deduplicate UMIs and generate the barcode sparse matrices for both.",
        "snakemake_target": "bcwithqc_count",
    },
    "count": {
        "help": "Create downstream count tables.",
        "description": "Create the downstream count tables and mapping summary outputs.",
        "snakemake_target": "count",
    },
    "MAUDE": {
        "help": "Run MAUDE analysis.",
        "description": "Run the MAUDE analysis on the downstream count data.",
        "snakemake_target": "MAUDE",
    },
    "plot": {
        "help": "Generate plots and summaries.",
        "description": "Generate the final summary plots and QC reports.",
        "snakemake_target": "plot",
    },
}


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

    # Give every run its own editable config.yaml.
    # Do not overwrite an existing config in case the user wants to edit it.
    config_src = source_root / "config.yaml"
    config_dst = pipeline_dir / "config.yaml"
    if not config_dst.exists():
        shutil.copy2(config_src, config_dst)

    # sample_profiles are always be copied, not symlinked,
    # so users can inspect or edit the local copy.
    sample_profiles_src = source_root / "sample_profiles"
    sample_profiles_dst = pipeline_dir / "sample_profiles"

    if not sample_profiles_src.exists():
        sys.exit(f"\nERROR\nMissing pipeline component:\n{sample_profiles_src}")

    if sample_profiles_dst.exists() or sample_profiles_dst.is_symlink():
        if sample_profiles_dst.is_symlink() or sample_profiles_dst.is_file():
            sample_profiles_dst.unlink()
        else:
            shutil.rmtree(sample_profiles_dst)

    shutil.copytree(sample_profiles_src, sample_profiles_dst)

    # Ensure all remaining pipeline components are present as symlinks.
    for item in PIPELINE_ITEMS:
        if item == "sample_profiles":
            continue

        src = source_root / item
        dest = pipeline_dir / item

        if not src.exists():
            sys.exit(f"\nERROR\nMissing pipeline component:\n{src}")

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
    if args.custom_profile:
        args.custom_profile = str(Path(args.custom_profile).resolve())


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
  
def set_nested(cfg, dotted_key, value):
    parts = dotted_key.split(".")
    cur = cfg
    for part in parts[:-1]:
        cur = cur.setdefault(part, {})
    cur[parts[-1]] = value


def apply_profile_overrides(profile_name, profile_config_path, args):
    PROFILE_RESOURCE_KEYS = {
        "local": {
            # No cluster-specific overrides are confirmed for local execution.
        },
        "slurm": {
            "account": "default-resources.slurm_account",
            "partition": "default-resources.slurm_partition",
            "qos": "default-resources.slurm_qos",
            "reservation": "default-resources.slurm_reservation",
        },
        "sge": {
            "partition": "default-resources.sge_queue",
            "account": "default-resources.sge_project",
            # "reservation": None,  # no documented SGE reservation setting
            # "qos": None,         # no documented SGE qos setting
            # "sge_pe" exists, but it is a separate parallel-environment setting,
            # not a generic account/partition/qos/reservation equivalent.
        },
        "lsf": {
            "partition": "default-resources.lsf_queue",
            "account": "default-resources.lsf_project",
            # "reservation": None,  # not confirmed in the docs
            # "qos": None,         # not confirmed in the docs
        },
        "htcondor": {
            # No confirmed account/partition/qos/reservation equivalents in the docs.
        },
        "cluster-generic": {
            # No confirmed account/partition/qos/reservation equivalents in the docs.
        },
    }
    with open(profile_config_path, "r") as handle:
        cfg = yaml.safe_load(handle) or {}

    mapping = PROFILE_RESOURCE_KEYS.get(profile_name, {})
    requested = {
        "account": args.account,
        "partition": args.partition,
        "qos": args.qos,
        "reservation": args.reservation,
    }

    unsupported = []

    for generic_name, value in requested.items():
        if value is None:
            continue

        if generic_name not in mapping:
            unsupported.append(generic_name)
            continue

        set_nested(cfg, mapping[generic_name], value)

    if unsupported:
        supported = ", ".join(sorted(mapping.keys())) if mapping else "none"
        raise ValueError(
            f"Profile '{profile_name}' does not support override(s): "
            + ", ".join(sorted(unsupported))
            + f". Supported overrides: {supported}"
        )

    with open(profile_config_path, "w") as handle:
        yaml.safe_dump(cfg, handle, sort_keys=False)
        

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

    return ",".join(f"{mount}:{mount}" for mount in mounts)


def add_common_arguments(parser):
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing FASTQ files.",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output master directory for the analysis.",
    )

    parser.add_argument(
        "--library-file",
        required=True,
        help="Library annotation file containing the guide sequences. See 'library construction' in the Documentation",
    )

    parser.add_argument(
        "--bcwithqc-config-file",
        required=True,
        help="bcwithqc JSON configuration file.",
    )

    parser.add_argument(
        "--experiment-info-file",
        required=True,
        help="Experiment information spreadsheet containing information for the fastq files and the bins. See 'Experiment info file Layout' in the Documentation",
    )

    parser.add_argument(
        "--reads-or-umis",
        choices=("reads", "umis"),
        default="reads",
        help="Use read counts or UMI counts (default: reads).",
    )
    
    parser.add_argument(
        "--profile",
        choices=("local", "cluster-generic", "slurm", "htcondor", "lsf", "sge"),
        default="local",
        help="Sets the environment/cluster that the snakemake pipeline is run on (default: local)",
    )
    
    parser.add_argument(
        "--custom-profile",
        default=None,
        help="Optional Snakemake profile path. If set, this overrides the built-in default profile.",
    )
    
    parser.add_argument("--account", default=None)
    parser.add_argument("--partition", default=None)
    parser.add_argument("--qos", default=None)
    parser.add_argument("--reservation", default=None)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="mecas",
        description="Run the MECAS amplicon barcode analysis pipeline.\n\nUse 'mecas COMMAND --help' for stage-specific help.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    subparsers = parser.add_subparsers(
        dest="command",
        # required=True, # This needs python 3.6
        metavar="COMMAND",
    )

    common_parser = argparse.ArgumentParser(add_help=False)
    add_common_arguments(common_parser)

    for command_name, command_spec in COMMAND_SPECS.items():
        command_parser = subparsers.add_parser(
            command_name,
            parents=[common_parser],
            help=command_spec["help"],
            description=command_spec["description"],
        )
        command_parser.set_defaults(snakemake_target=command_spec["snakemake_target"])

    version_parser = subparsers.add_parser(
        "version",
        help="Print MECAS version information and exit.",
        description="Print MECAS version information and exit.",
    )
    version_parser.set_defaults(version_command=True)

    return parser


def main():
    parser = build_parser()

    # Forward all unknown options directly to Snakemake.
    args, snakemake_args = parser.parse_known_args()

    if getattr(args, "version_command", False):
        print(f"mecas {VERSION}")
        raise SystemExit(0)
    
    if not getattr(args, "command", None):
        parser.print_help()
        raise SystemExit(2)

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
    
    # Use either a custom provided profile or one of the default ones
    default_profile = Path("sample_profiles") / args.profile
    profile_arg = args.custom_profile if args.custom_profile else default_profile
    
    # Aplly profile overrides
    if not args.custom_profile:
        apply_profile_overrides(
            args.profile,
            pipeline_dir / "sample_profiles" / args.profile / "config.yaml",
            args,
        )
        
    cmd = [
        "snakemake",
        args.snakemake_target,
        "--sdm",
        "apptainer",
        "--apptainer-args",
        f"--bind {apptainer_bind_string}",
        "--profile",
        str(profile_arg),
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
