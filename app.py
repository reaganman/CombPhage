#!/usr/bin/env python3
from __future__ import annotations

import re
import shutil
from datetime import datetime
from pathlib import Path

import pandas as pd
import streamlit as st

from combphage_pipeline import CombPhagePipeline, CombPhagePipelineError


REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_RUNS_DIR = REPO_ROOT / "combphage_runs"

STAGES = ["fetch", "prepare", "overlaps", "filter", "visualize"]
STAGE_LABELS = {
    "fetch": "GenBank input",
    "prepare": "Prepare genomes",
    "overlaps": "Cluster proteins and find overlaps",
    "filter": "Filter overlaps",
    "visualize": "Generate visualization",
}


def safe_name(value):
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    value = value.strip("._")
    return value or "combphage_run"


def parse_accessions(text):
    return [value for value in re.split(r"[\s,;]+", text.strip()) if value]


def save_uploads(files, destination):
    destination.mkdir(parents=True, exist_ok=True)
    used = set()
    for uploaded in files:
        filename = Path(uploaded.name).name
        stem = Path(filename).stem
        suffix = Path(filename).suffix
        candidate = filename
        counter = 2
        while candidate in used:
            candidate = f"{stem}_{counter}{suffix}"
            counter += 1
        used.add(candidate)
        (destination / candidate).write_bytes(uploaded.getvalue())


def read_tsv(path):
    if path.exists() and path.stat().st_size > 0:
        return pd.read_csv(path, sep="\t")
    return None


def show_results(run_dir):
    overlaps_dir = run_dir / "feature_overlaps"

    filtered = read_tsv(overlaps_dir / "filtered_overlaps.tsv")
    membership = read_tsv(overlaps_dir / "cluster_membership.tsv")
    resolved = read_tsv(overlaps_dir / "resolved_clusters.tsv")
    ambiguous = read_tsv(overlaps_dir / "ambiguous_clusters.tsv")

    st.subheader("Results")
    cols = st.columns(4)
    cols[0].metric(
        "Final clusters",
        membership["cluster"].nunique()
        if membership is not None and "cluster" in membership else 0,
    )
    cols[1].metric("Filtered overlaps", len(filtered) if filtered is not None else 0)
    cols[2].metric("Resolved CDS assignments", len(resolved) if resolved is not None else 0)
    cols[3].metric("Raw ambiguous CDSs", len(ambiguous) if ambiguous is not None else 0)

    table_specs = [
        ("Filtered overlaps", "filtered_overlaps.tsv"),
        ("Filter report", "overlap_filter_report.tsv"),
        ("Cluster membership", "cluster_membership.tsv"),
        ("Resolved clusters", "resolved_clusters.tsv"),
        ("Ambiguous clusters", "ambiguous_clusters.tsv"),
        ("Incomplete clusters", "incomplete_clusters.tsv"),
        ("Unresolved members", "unresolved_ambiguous_members.tsv"),
        ("Unclustered CDSs", "unclustered_cds.tsv"),
    ]

    available = [
        (label, overlaps_dir / filename)
        for label, filename in table_specs
        if (overlaps_dir / filename).exists()
    ]

    if available:
        tabs = st.tabs([label for label, _ in available])
        for tab, (label, path) in zip(tabs, available):
            with tab:
                df = read_tsv(path)
                if df is not None:
                    st.dataframe(df, use_container_width=True, hide_index=True)
                    st.download_button(
                        f"Download {path.name}",
                        data=path.read_bytes(),
                        file_name=path.name,
                        mime="text/tab-separated-values",
                        key=f"download_{path.name}",
                    )

    st.subheader("Genome visualization")
    viz_dir = overlaps_dir / "lovis4u_overlaps"
    if viz_dir.exists():
        raster_files = sorted(
            list(viz_dir.glob("*.png"))
            + list(viz_dir.glob("*.jpg"))
            + list(viz_dir.glob("*.jpeg"))
        )
        for path in raster_files:
            st.image(str(path), caption=path.name)

        export_files = sorted(
            p for p in viz_dir.iterdir()
            if p.is_file()
            and p.suffix.lower() in {".pdf", ".svg", ".png", ".jpg", ".jpeg"}
        )
        for path in export_files:
            st.download_button(
                f"Download {path.name}",
                data=path.read_bytes(),
                file_name=path.name,
                key=f"viz_{path.name}",
            )
    else:
        st.info("No LoVis4u visualization was generated.")

    archive_base = run_dir.parent / f"{run_dir.name}_results"
    archive_path = Path(shutil.make_archive(str(archive_base), "zip", root_dir=run_dir))
    st.download_button(
        "Download complete results (.zip)",
        data=archive_path.read_bytes(),
        file_name=archive_path.name,
        mime="application/zip",
    )

    with st.expander("Pipeline logs"):
        logs_dir = run_dir / "logs"
        for path in sorted(logs_dir.glob("*.log")):
            st.markdown(f"**{path.name}**")
            st.code(path.read_text(encoding="utf-8", errors="replace")[-12000:])


st.set_page_config(page_title="CombPhage", layout="wide")
st.title("CombPhage")
st.caption("Identify conserved feature boundary overlaps for designing interchangeable phage genome fragments.")

with st.sidebar:
    st.header("Analysis settings")
    min_overlap = st.number_input(
        "Minimum overlap size (bp)", min_value=1, max_value=500, value=25, step=1
    )
    max_homopolymer = st.number_input(
        "Maximum allowed homopolymer length", min_value=1, max_value=20, value=4, step=1
    )
    generate_viz = st.checkbox("Generate LoVis4u visualization", value=True)
    runs_dir_text = st.text_input("Runs directory", value=str(DEFAULT_RUNS_DIR))

input_mode = st.radio(
    "Genome input",
    ["Upload GenBank files", "Fetch from NCBI accessions"],
    horizontal=True,
)

uploaded_files = []
accession_text = ""
email = ""

if input_mode == "Upload GenBank files":
    uploaded_files = st.file_uploader(
        "Upload annotated GenBank files",
        type=["gb", "gbk", "gbff", "genbank"],
        accept_multiple_files=True,
    )
else:
    accession_text = st.text_area(
        "NCBI nucleotide accessions",
        placeholder="MZ501081.1 MZ501078.1 V01146.1",
    )
    email = st.text_input("Email for NCBI Entrez", placeholder="email@example.com")

run_name = st.text_input(
    "Run name",
    value="combphage_" + datetime.now().strftime("%Y%m%d_%H%M%S"),
)

if st.button("Run CombPhage", type="primary", use_container_width=True):
    error = None

    if input_mode == "Upload GenBank files":
        if len(uploaded_files) < 2:
            error = "Upload at least two GenBank files."
    else:
        accessions = parse_accessions(accession_text)
        if len(accessions) < 2:
            error = "Provide at least two NCBI accessions."
        elif not email.strip():
            error = "Provide an email address for NCBI Entrez."

    if error:
        st.error(error)
    else:
        runs_dir = Path(runs_dir_text).expanduser().resolve()
        run_dir = runs_dir / safe_name(run_name)

        if run_dir.exists():
            st.error(f"Run directory already exists: {run_dir}\nChoose a different run name.")
        else:
            stage_states = {stage: "pending" for stage in STAGES}
            progress = st.progress(0)
            status = st.empty()
            stage_display = st.empty()

            def redraw():
                icons = {
                    "pending": "[ ]",
                    "running": "[>]",
                    "complete": "[x]",
                    "skipped": "[-]",
                    "error": "[!]",
                }
                stage_display.code(
                    "\n".join(
                        f"{icons[stage_states[s]]} {STAGE_LABELS[s]}"
                        for s in STAGES
                    )
                )

            def callback(stage, state, message):
                stage_states[stage] = state
                status.write(message)
                redraw()
                done = sum(
                    s in {"complete", "skipped"}
                    for s in stage_states.values()
                )
                progress.progress(done / len(STAGES))

            redraw()

            try:
                pipeline = CombPhagePipeline(REPO_ROOT, run_dir, callback)

                if input_mode == "Upload GenBank files":
                    incoming = run_dir / "_uploaded_genbanks"
                    save_uploads(uploaded_files, incoming)
                    pipeline.use_local_genbanks(incoming)
                else:
                    pipeline.fetch_genbanks(parse_accessions(accession_text), email)

                pipeline.run_all(
                    min_size=int(min_overlap),
                    max_homo_size=int(max_homopolymer),
                    visualize=generate_viz,
                )

                progress.progress(1.0)
                status.success("CombPhage analysis complete.")
                st.session_state["last_combphage_run"] = str(run_dir)

            except CombPhagePipelineError as exc:
                st.error(str(exc))
                if exc.log_file and exc.log_file.exists():
                    with st.expander(f"{exc.stage} log", expanded=True):
                        st.code(
                            exc.log_file.read_text(
                                encoding="utf-8", errors="replace"
                            )[-12000:]
                        )
            except Exception as exc:
                st.exception(exc)

last_run = st.session_state.get("last_combphage_run")
if last_run:
    run_path = Path(last_run)
    if run_path.exists():
        st.divider()
        show_results(run_path)
