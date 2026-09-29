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
