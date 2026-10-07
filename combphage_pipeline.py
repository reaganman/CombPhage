#!/usr/bin/env python3
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


class CombPhagePipelineError(RuntimeError):
    def __init__(self, stage, message, log_file=None):
        self.stage = stage
        self.log_file = log_file
        super().__init__(message)


class CombPhagePipeline:
    def __init__(self, repo_root, run_dir, progress_callback=None):
        self.repo_root = Path(repo_root).resolve()
        self.scripts_dir = self.repo_root / "Scripts"
        self.run_dir = Path(run_dir).resolve()
        self.progress_callback = progress_callback

        self.genbanks_dir = self.run_dir / "genbanks"
        self.prepped_dir = self.run_dir / "prepped_genbanks"
        self.overlaps_dir = self.run_dir / "feature_overlaps"
        self.logs_dir = self.run_dir / "logs"

        self.run_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        self._validate_scripts()

    def _validate_scripts(self):
        required = [
            "fetch_genbanks.py",
            "prepare_genbanks.py",
            "find_overlaps.py",
            "filter_overlaps.py",
            "visualize.py",
        ]
        missing = [name for name in required if not (self.scripts_dir / name).exists()]
        if missing:
            raise FileNotFoundError(
                "Missing required CombPhage script(s): "
                + ", ".join(missing)
                + f"\nExpected them under: {self.scripts_dir}"
            )

    def _notify(self, stage, state, message):
        if self.progress_callback:
            self.progress_callback(stage, state, message)

    def _run_command(self, stage, cmd):
        log_file = self.logs_dir / f"{stage}.log"
        self._notify(stage, "running", f"Running {stage}...")

        with log_file.open("w", encoding="utf-8") as log:
            log.write("$ " + " ".join(str(x) for x in cmd) + "\n\n")
            log.flush()
            process = subprocess.run(
                [str(x) for x in cmd],
                cwd=self.repo_root,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
            )

        if process.returncode != 0:
            self._notify(stage, "error", f"{stage} failed.")
            raise CombPhagePipelineError(
                stage,
                f"CombPhage stage '{stage}' failed with exit code {process.returncode}.",
                log_file,
            )

        self._notify(stage, "complete", f"{stage} complete.")

    def fetch_genbanks(self, accessions, email):
        accessions = [str(x).strip() for x in accessions if str(x).strip()]
        if not accessions:
            raise ValueError("At least one accession is required.")
        if not email.strip():
            raise ValueError("An email address is required for NCBI fetching.")

        self.genbanks_dir.mkdir(parents=True, exist_ok=True)
        self._run_command(
            "fetch",
            [
                sys.executable,
                self.scripts_dir / "fetch_genbanks.py",
                "-a",
                *accessions,
                "-o",
                self.genbanks_dir,
                "--email",
                email.strip(),
            ],
        )

    def use_local_genbanks(self, source_dir):
        source_dir = Path(source_dir).resolve()
        extensions = {".gb", ".gbk", ".gbff", ".genbank"}
        files = [
            p for p in source_dir.iterdir()
            if p.is_file() and p.suffix.lower() in extensions
        ]
        if not files:
            raise ValueError(f"No GenBank files found in {source_dir}")

        self.genbanks_dir.mkdir(parents=True, exist_ok=True)
        for src in files:
            dest = self.genbanks_dir / src.name
            if src.resolve() != dest.resolve():
                shutil.copy2(src, dest)

        self._notify("fetch", "skipped", f"Using {len(files)} uploaded GenBank file(s).")

    def prepare_genbanks(self):
        self._run_command(
            "prepare",
            [
                sys.executable,
                self.scripts_dir / "prepare_genbanks.py",
                "-i",
                self.genbanks_dir,
                "-o",
                self.prepped_dir,
            ],
        )

    def find_overlaps(self):
        self._run_command(
            "overlaps",
            [
                sys.executable,
                self.scripts_dir / "find_overlaps.py",
                "--input-dir",
                self.prepped_dir,
                "--outdir",
                self.overlaps_dir,
            ],
        )

    def filter_overlaps(self, min_size=25, max_homo_size=4):
        self._run_command(
            "filter",
            [
                sys.executable,
                self.scripts_dir / "filter_overlaps.py",
                "--overlaps",
                self.overlaps_dir / "cluster_overlap_summary.tsv",
                "--outdir",
                self.overlaps_dir,
                "--min_size",
                str(min_size),
                "--max_homo_size",
                str(max_homo_size),
            ],
        )

    def visualize(self):
        self._run_command(
            "visualize",
            [
                sys.executable,
                self.scripts_dir / "visualize.py",
                "--genbanks",
                self.prepped_dir / "genbanks",
                "--overlaps-dir",
                self.overlaps_dir,
                "--force",
            ],
        )

    def run_all(self, min_size=25, max_homo_size=4, visualize=True):
        self.prepare_genbanks()
        self.find_overlaps()
        self.filter_overlaps(min_size, max_homo_size)
        if visualize:
            self.visualize()
        else:
            self._notify("visualize", "skipped", "Visualization skipped.")


def _cli_progress(stage, state, message):
    """Simple terminal progress reporter used by the command-line interface."""
    symbols = {
        "running": "[RUN]",
        "complete": "[OK]",
        "skipped": "[SKIP]",
        "error": "[ERROR]",
    }
    print(f"{symbols.get(state, '[INFO]')} {message}", flush=True)


def build_parser():
    """Build the command-line parser without affecting GUI imports."""
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Run the complete CombPhage overlap-discovery pipeline from "
            "GenBank files or NCBI nucleotide accessions."
        )
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--genbanks",
        help="Directory containing annotated GenBank files.",
    )
    input_group.add_argument(
        "--accessions",
        nargs="+",
        help=(
            "NCBI nucleotide accession(s) to fetch before analysis, e.g. "
            "MZ501081.1 MZ501078.1 V01146.1"
        ),
    )

    parser.add_argument(
        "--email",
        help="Email address required by NCBI Entrez when using --accessions.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory for the complete CombPhage run.",
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        default=25,
        help="Minimum filtered overlap size in bp (default: 25).",
    )
    parser.add_argument(
        "--max-homopolymer",
        type=int,
        default=4,
        help="Maximum allowed homopolymer length (default: 4).",
    )
    parser.add_argument(
        "--no-visualization",
        action="store_true",
        help="Skip the final LoVis4u visualization stage.",
    )
    parser.add_argument(
        "--repo-root",
        default=str(Path(__file__).resolve().parent),
        help=(
            "CombPhage repository root containing the Scripts directory. "
            "Defaults to the directory containing this file."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Remove an existing non-empty output directory before starting. "
            "Use with care."
        ),
    )

    return parser


def main():
    """Command-line entry point. The class above remains importable by app.py."""
    parser = build_parser()
    args = parser.parse_args()

    if args.accessions:
        if len(args.accessions) < 2:
            parser.error("CombPhage requires at least two genomes for comparison.")
        if not args.email:
            parser.error("--email is required when using --accessions.")

    if args.min_overlap < 1:
        parser.error("--min-overlap must be at least 1.")

    if args.max_homopolymer < 1:
        parser.error("--max-homopolymer must be at least 1.")

    run_dir = Path(args.outdir).expanduser().resolve()

    if run_dir.exists() and any(run_dir.iterdir()):
        if args.force:
            shutil.rmtree(run_dir)
        else:
            parser.error(
                f"Output directory already exists and is not empty: {run_dir}\n"
                "Choose a new directory or rerun with --force."
            )

    pipeline = CombPhagePipeline(
        repo_root=args.repo_root,
        run_dir=run_dir,
        progress_callback=_cli_progress,
    )

    try:
        if args.accessions:
            _cli_progress(
                "fetch",
                "running",
                f"Fetching {len(args.accessions)} GenBank record(s) from NCBI...",
            )
            pipeline.fetch_genbanks(args.accessions, args.email)
        else:
            source_dir = Path(args.genbanks).expanduser().resolve()
            if not source_dir.exists() or not source_dir.is_dir():
                parser.error(f"GenBank directory does not exist: {source_dir}")

            pipeline.use_local_genbanks(source_dir)

        pipeline.run_all(
            min_size=args.min_overlap,
            max_homo_size=args.max_homopolymer,
            visualize=not args.no_visualization,
        )

    except CombPhagePipelineError as exc:
        print(f"\nCombPhage failed during stage: {exc.stage}", file=sys.stderr)
        print(str(exc), file=sys.stderr)
        if exc.log_file:
            print(f"Log: {exc.log_file}", file=sys.stderr)
        return 1
    except Exception as exc:
        print(f"\nCombPhage failed: {exc}", file=sys.stderr)
        return 1

    print("\nCombPhage analysis complete.")
    print(f"Results: {run_dir}")
    print(f"Prepared genomes: {pipeline.prepped_dir}")
    print(f"Overlap results: {pipeline.overlaps_dir}")
    print(f"Logs: {pipeline.logs_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
