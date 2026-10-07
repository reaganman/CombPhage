#!/usr/bin/env python3
from __future__ import annotations

import json
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import streamlit as st

from combphage_pipeline import CombPhagePipeline, CombPhagePipelineError


REPO_ROOT = Path(__file__).resolve().parent
SCRIPTS_DIR = REPO_ROOT / "Scripts"
DEFAULT_RUNS_DIR = REPO_ROOT / "combphage_runs"

# design_assemblies.py lives with the other CombPhage pipeline scripts.
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from design_assemblies import (  # noqa: E402
    DEFAULT_MAX_CHIMERAS,
    DEFAULT_MAX_COMBINATIONS,
    DEFAULT_MAX_FRAGMENT_SIZE,
    DEFAULT_IDEAL_OVERLAP_LENGTH,
    DEFAULT_IDEAL_GC,
    DEFAULT_TOP_N,
    DEFAULT_WEIGHT_GROUP_STD_DEV,
    DEFAULT_WEIGHT_INDIVIDUAL_PENALTY,
    DEFAULT_SHORT_OVERLAP_MULTIPLIER,
    DEFAULT_LONG_OVERLAP_MULTIPLIER,
    DEFAULT_GC_MULTIPLIER,
    DEFAULT_LARGE_FRAGMENT_MULTIPLIER,
    candidate_summary_table,
    design_from_candidates,
    materialize_design,
    prepare_design_context,
    rank_designs,
    candidate_from_subsequence,
    write_design_outputs,
    write_single_design,
)


from analyze_design import (  # noqa: E402
    DEFAULT_MAX_PAIRWISE_ASSEMBLIES,
    DEFAULT_TOP_VARIABLE_CDS,
    analyze_fragment_design,
    analyze_materialized_design,
)


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


def archive_directory(directory: Path, archive_name: str) -> Path:
    archive_base = directory.parent / archive_name
    archive_path = Path(
        shutil.make_archive(
            str(archive_base),
            "zip",
            root_dir=directory,
        )
    )
    return archive_path


def show_results(run_dir):
    overlaps_dir = run_dir / "feature_overlaps"

    filtered = read_tsv(overlaps_dir / "filtered_overlaps.tsv")
    membership = read_tsv(overlaps_dir / "cluster_membership.tsv")
    resolved = read_tsv(overlaps_dir / "resolved_clusters.tsv")
    ambiguous = read_tsv(overlaps_dir / "ambiguous_clusters.tsv")

    st.subheader("Overlap analysis results")

    cols = st.columns(4)
    cols[0].metric(
        "Final clusters",
        membership["cluster"].nunique()
        if membership is not None and "cluster" in membership
        else 0,
    )
    cols[1].metric(
        "Filtered overlaps",
        len(filtered) if filtered is not None else 0,
    )
    cols[2].metric(
        "Resolved CDS assignments",
        len(resolved) if resolved is not None else 0,
    )
    cols[3].metric(
        "Raw ambiguous CDSs",
        len(ambiguous) if ambiguous is not None else 0,
    )

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

    with st.expander("Overlap result tables", expanded=False):
        if available:
            tabs = st.tabs([label for label, _ in available])

            for tab, (label, path) in zip(tabs, available):
                with tab:
                    df = read_tsv(path)
                    if df is not None:
                        st.dataframe(
                            df,
                            use_container_width=True,
                            hide_index=True,
                        )
                        st.download_button(
                            f"Download {path.name}",
                            data=path.read_bytes(),
                            file_name=path.name,
                            mime="text/tab-separated-values",
                            key=f"download_{run_dir.name}_{path.name}",
                        )
        else:
            st.info("No overlap result tables were found.")

    with st.expander("Overlap visualization", expanded=False):
        viz_dir = overlaps_dir / "lovis4u_overlaps"

        if not viz_dir.exists():
            st.info("No LoVis4u visualization was generated.")
        else:
            raster_files = sorted(
                list(viz_dir.glob("*.png"))
                + list(viz_dir.glob("*.jpg"))
                + list(viz_dir.glob("*.jpeg"))
            )
            pdf_files = sorted(viz_dir.glob("*.pdf"))

            if raster_files:
                for path in raster_files:
                    st.image(str(path), caption=path.name)
            elif pdf_files:
                # Streamlit 1.58+ can render PDFs directly when the optional
                # streamlit-pdf dependency is installed.
                try:
                    st.pdf(str(pdf_files[0]), height=750)
                except Exception:
                    st.info(
                        "The visualization PDF was generated, but inline PDF "
                        "preview is unavailable in this Streamlit environment. "
                        "Install the PDF extra with `pip install \"streamlit[pdf]\"` "
                        "to enable it."
                    )
            else:
                st.info("The LoVis4u output directory contains no previewable image or PDF.")

            export_files = sorted(
                path
                for path in viz_dir.iterdir()
                if path.is_file()
                and path.suffix.lower()
                in {".pdf", ".svg", ".png", ".jpg", ".jpeg"}
            )

            for path in export_files:
                st.download_button(
                    f"Download {path.name}",
                    data=path.read_bytes(),
                    file_name=path.name,
                    key=f"viz_{run_dir.name}_{path.name}",
                )

    with st.expander("Pipeline logs", expanded=False):
        logs_dir = run_dir / "logs"
        if logs_dir.exists():
            for path in sorted(logs_dir.glob("*.log")):
                st.markdown(f"**{path.name}**")
                st.code(
                    path.read_text(
                        encoding="utf-8",
                        errors="replace",
                    )[-12000:]
                )

def assembly_prefix(run_dir: Path) -> str:
    return "assembly_" + safe_name(run_dir.name)


def reset_working_design(prefix: str, design: dict, primary_genome: str):
    ordered = sorted(
        design["combo"],
        key=lambda overlap: overlap["coords"][primary_genome]["center"],
    )

    selected_ids = [overlap["id"] for overlap in ordered]
    st.session_state[f"{prefix}_selected_ids"] = selected_ids

    for index, overlap in enumerate(ordered):
        candidate_id = overlap["id"]
        st.session_state[f"{prefix}_junction_{index}"] = candidate_id
        st.session_state[f"{prefix}_junction_{index}_previous"] = candidate_id
        # Rank 1 already contains the automatically optimized subsequence.
        # Keep that as the default working sequence while exposing the full
        # filtered overlap for arbitrary manual subsequence selection.
        st.session_state[f"{prefix}_junction_{index}_subsequence"] = overlap["seq"]

def ranked_summary_dataframe(designs, genomes):
    rows = []

    for design in designs:
        row = {
            "rank": design["rank"],
            "score": round(design["score"], 2),
            "junctions": "; ".join(
                overlap["id"] for overlap in design["combo"]
            ),
            "mean fragment-size SD": round(
                design["group_std_dev"],
                1,
            ),
            "overlap penalty": round(
                design["average_overlap_penalty"],
                2,
            ),
            "large-fragment penalty": round(
                design["large_fragment_penalty"],
                1,
            ),
        }

        for genome in genomes:
            row[f"{genome} fragment sizes"] = "; ".join(
                f"{size:.0f}"
                for size in design["fragment_sizes"][genome]
            )

        rows.append(row)

    return pd.DataFrame(rows)


def current_design_tables(design, context):
    primary = context["primary_genome"]

    ordered = sorted(
        design["combo"],
        key=lambda overlap: overlap["coords"][primary]["center"],
    )

    junction_rows = []

    for order, overlap in enumerate(ordered, start=1):
        junction_rows.append(
            {
                "order": order,
                "candidate": overlap["id"],
                "cluster": overlap["cluster"],
                "region": overlap["region"],
                f"position in {primary}": round(
                    overlap["coords"][primary]["center"],
                    1,
                ),
                "length": overlap["length"],
                "GC %": round(overlap["gc"] * 100, 1),
                "working sequence": overlap["seq"],
            }
        )

    fragment_rows = []
    fragment_count = max(
        len(sizes)
        for sizes in design["fragment_sizes"].values()
    )

    for genome in context["genomes"]:
        row = {"genome": genome}
        sizes = design["fragment_sizes"][genome]

        for index in range(fragment_count):
            row[f"Fragment {index + 1}"] = (
                round(sizes[index]) if index < len(sizes) else None
            )

        row["size SD"] = round(design["std_devs"][genome], 1)
        fragment_rows.append(row)

    return pd.DataFrame(junction_rows), pd.DataFrame(fragment_rows)


def _read_analysis_overview(analysis_dir: Path):
    path = analysis_dir / "analysis_summary.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def _show_analysis_plots(analysis_dir: Path, key_prefix: str):
    plots_dir = analysis_dir / "plots"
    if not plots_dir.exists():
        return

    plot_files = sorted(plots_dir.glob("*.png"))
    if not plot_files:
        return

    st.markdown("**Visualizations**")
    for path in plot_files:
        st.image(str(path), caption=path.stem.replace("_", " "))


def show_fragment_analysis_outputs(run_dir: Path, prefix: str):
    working_dir = run_dir / "assembly_designs" / "working_design"
    analysis_dir = working_dir / "analysis"
    if not analysis_dir.exists():
        return

    overview = _read_analysis_overview(analysis_dir) or {}

    with st.expander("Current fragment analysis", expanded=True):
        metrics = st.columns(4)
        metrics[0].metric(
            "Fragments",
            overview.get("n_fragments", 0),
        )
        mean_frag = overview.get("mean_fragment_pairwise_nucleotide_identity")
        metrics[1].metric(
            "Mean fragment nt identity",
            f"{mean_frag:.2f}%" if mean_frag is not None else "—",
        )
        mean_cds = overview.get("mean_cds_cluster_nucleotide_identity")
        metrics[2].metric(
            "Mean CDS nt identity",
            f"{mean_cds:.2f}%" if mean_cds is not None else "—",
        )
        mean_prot = overview.get("mean_cds_cluster_protein_identity")
        metrics[3].metric(
            "Mean protein identity",
            f"{mean_prot:.2f}%" if mean_prot is not None else "—",
        )

        variable = overview.get("most_variable_fragment")
        if variable:
            st.info(
                f"Most variable fragment: Fragment "
                f"{variable['fragment_number']} "
                f"({variable['mean_pairwise_nucleotide_identity']:.2f}% "
                "mean pairwise nucleotide identity)."
            )

        table_specs = [
            ("Fragment summary", "fragment_summary.tsv"),
            ("Fragment pairwise identity", "fragment_pairwise_identity.tsv"),
            ("CDS cluster summary", "cds_cluster_summary.tsv"),
            ("CDS pairwise identity", "cds_pairwise_identity.tsv"),
        ]
        available = [
            (label, analysis_dir / filename)
            for label, filename in table_specs
            if (analysis_dir / filename).exists()
        ]
        if available:
            tabs = st.tabs([label for label, _ in available])
            for tab, (label, path) in zip(tabs, available):
                with tab:
                    df = read_tsv(path)
                    if df is not None:
                        st.dataframe(df, use_container_width=True, hide_index=True)

        _show_analysis_plots(analysis_dir, f"{prefix}_fragment")

        archive_path = archive_directory(
            working_dir,
            "working_design_fragment_analysis",
        )
        st.download_button(
            "Download current fragment analysis (.zip)",
            data=archive_path.read_bytes(),
            file_name=archive_path.name,
            mime="application/zip",
            key=f"{prefix}_fragment_analysis_download",
        )


def show_full_design_analysis_outputs(run_dir: Path, prefix: str):
    selected_dir = run_dir / "assembly_designs" / "selected_design"
    analysis_dir = selected_dir / "analysis"
    if not analysis_dir.exists():
        return

    overview = _read_analysis_overview(analysis_dir) or {}

    with st.expander("Full assembly-library analysis", expanded=True):
        metrics = st.columns(4)
        metrics[0].metric("Chimeras", overview.get("n_chimeras", 0))
        mean_frag = overview.get("mean_fragment_pairwise_nucleotide_identity")
        metrics[1].metric(
            "Mean fragment nt identity",
            f"{mean_frag:.2f}%" if mean_frag is not None else "—",
        )
        mean_cds = overview.get("mean_cds_cluster_nucleotide_identity")
        metrics[2].metric(
            "Mean CDS nt identity",
            f"{mean_cds:.2f}%" if mean_cds is not None else "—",
        )
        nearest = overview.get("mean_chimera_nearest_parent_identity")
        metrics[3].metric(
            "Mean chimera nearest-parent identity",
            f"{nearest:.2f}%" if nearest is not None else "—",
        )

        table_specs = [
            ("Fragment summary", "fragment_summary.tsv"),
            ("CDS cluster summary", "cds_cluster_summary.tsv"),
            ("Assembly summary", "assembly_summary.tsv"),
            ("Assembly vs parent", "assembly_vs_parent_identity.tsv"),
            ("Assembly pairwise identity", "assembly_pairwise_identity.tsv"),
        ]
        available = [
            (label, analysis_dir / filename)
            for label, filename in table_specs
            if (analysis_dir / filename).exists()
        ]
        if available:
            tabs = st.tabs([label for label, _ in available])
            for tab, (label, path) in zip(tabs, available):
                with tab:
                    df = read_tsv(path)
                    if df is not None:
                        st.dataframe(df, use_container_width=True, hide_index=True)

        _show_analysis_plots(analysis_dir, f"{prefix}_full")

        chimera_log = read_tsv(
            selected_dir / "chimeras" / "chimera_assembly_log.tsv"
        )
        if chimera_log is not None:
            with st.expander("Chimera assembly log", expanded=False):
                st.dataframe(
                    chimera_log,
                    use_container_width=True,
                    hide_index=True,
                )

        archive_path = archive_directory(
            selected_dir,
            "selected_design_full_analysis",
        )
        st.download_button(
            "Download FASTAs and full analysis (.zip)",
            data=archive_path.read_bytes(),
            file_name=archive_path.name,
            mime="application/zip",
            key=f"{prefix}_full_analysis_download",
        )


def initialize_assembly_design(run_dir: Path, settings: dict):
    """Rank assembly designs and load rank 1 as the editable working design."""
    overlaps_file = run_dir / "feature_overlaps" / "filtered_overlaps.tsv"
    fasta_file = run_dir / "prepped_genbanks" / "genomes.fasta"
    prefix = assembly_prefix(run_dir)
    outdir = run_dir / "assembly_designs"

    context = prepare_design_context(
        overlaps=overlaps_file,
        fasta=fasta_file,
        ideal_overlap_length=int(settings["ideal_overlap_length"]),
        ideal_gc=float(settings["ideal_gc_percent"]) / 100.0,
        short_overlap_multiplier=float(settings["short_overlap_multiplier"]),
        long_overlap_multiplier=float(settings["long_overlap_multiplier"]),
        gc_multiplier=float(settings["gc_multiplier"]),
    )

    designs, total_combinations = rank_designs(
        context,
        num_fragments=int(settings["num_fragments"]),
        topology=settings["topology"],
        top_n=int(settings["top_n"]),
        max_fragment_size=int(settings["max_fragment_size"]),
        max_combinations=int(settings["max_combinations"]),
        weight_group_std_dev=float(settings["weight_group_std_dev"]),
        weight_individual_penalty=float(settings["weight_individual_penalty"]),
        large_fragment_multiplier=float(settings["large_fragment_multiplier"]),
    )

    write_design_outputs(
        designs,
        context["genomes"],
        outdir / "ranked_designs",
    )

    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "assembly_design_settings.json").write_text(
        json.dumps(settings, indent=2)
    )

    st.session_state[f"{prefix}_context"] = context
    st.session_state[f"{prefix}_designs"] = designs
    st.session_state[f"{prefix}_total_combinations"] = total_combinations
    st.session_state[f"{prefix}_settings"] = settings
    st.session_state[f"{prefix}_generated_signature"] = None
    st.session_state[f"{prefix}_fragment_analysis_signature"] = None
    st.session_state[f"{prefix}_full_analysis_signature"] = None

    reset_working_design(
        prefix,
        designs[0],
        context["primary_genome"],
    )

    return context, designs, total_combinations


def ensure_assembly_design_state(run_dir: Path):
    """Rehydrate assembly-design state after a Streamlit/server restart."""
    prefix = assembly_prefix(run_dir)
    context = st.session_state.get(f"{prefix}_context")
    designs = st.session_state.get(f"{prefix}_designs")
    settings = st.session_state.get(f"{prefix}_settings")

    if context and designs and settings:
        return context, designs, settings

    settings_file = run_dir / "assembly_designs" / "assembly_design_settings.json"
    if not settings_file.exists():
        return None, None, None

    try:
        settings = json.loads(settings_file.read_text())
        context, designs, _ = initialize_assembly_design(run_dir, settings)
        return context, designs, settings
    except Exception as exc:
        st.warning(f"Could not reload assembly-design state: {exc}")
        return None, None, None


def _default_design_settings(run_dir: Path):
    defaults = {
        "num_fragments": 4,
        "topology": "circular",
        "top_n": DEFAULT_TOP_N,
        "ideal_overlap_length": DEFAULT_IDEAL_OVERLAP_LENGTH,
        "ideal_gc_percent": DEFAULT_IDEAL_GC * 100.0,
        "max_fragment_size": DEFAULT_MAX_FRAGMENT_SIZE,
        "min_overlap": 1,
        "weight_group_std_dev": DEFAULT_WEIGHT_GROUP_STD_DEV,
        "weight_individual_penalty": DEFAULT_WEIGHT_INDIVIDUAL_PENALTY,
        "short_overlap_multiplier": DEFAULT_SHORT_OVERLAP_MULTIPLIER,
        "long_overlap_multiplier": DEFAULT_LONG_OVERLAP_MULTIPLIER,
        "gc_multiplier": DEFAULT_GC_MULTIPLIER,
        "large_fragment_multiplier": DEFAULT_LARGE_FRAGMENT_MULTIPLIER,
        "max_combinations": DEFAULT_MAX_COMBINATIONS,
    }

    settings_file = (
        run_dir / "assembly_designs" / "assembly_design_settings.json"
    )
    if settings_file.exists():
        try:
            defaults.update(json.loads(settings_file.read_text()))
        except Exception:
            pass

    return defaults


def show_assembly_design(run_dir: Path):
    overlaps_file = run_dir / "feature_overlaps" / "filtered_overlaps.tsv"
    fasta_file = run_dir / "prepped_genbanks" / "genomes.fasta"

    if not overlaps_file.exists() or not fasta_file.exists():
        return

    prefix = assembly_prefix(run_dir)
    outdir = run_dir / "assembly_designs"
    defaults = _default_design_settings(run_dir)

    st.divider()
    st.header("Assembly design")
    st.write(
        "Rank candidate overlap sets using the completed CombPhage overlap "
        "analysis. These settings can be changed and re-ranked without rerunning "
        "genome preparation, protein clustering, overlap discovery, filtering, "
        "or visualization."
    )

    with st.expander("Assembly design settings", expanded=True):
        col1, col2, col3 = st.columns(3)
        num_fragments = col1.number_input(
            "Number of fragments",
            min_value=2,
            max_value=20,
            value=int(defaults["num_fragments"]),
            step=1,
            key=f"{prefix}_setting_num_fragments",
        )
        topology_label = col2.selectbox(
            "Genome topology",
            ["Circular", "Linear"],
            index=0 if defaults["topology"] == "circular" else 1,
            key=f"{prefix}_setting_topology",
        )
        top_n = col3.number_input(
            "Ranked designs to retain",
            min_value=1,
            max_value=20,
            value=int(defaults["top_n"]),
            step=1,
            key=f"{prefix}_setting_top_n",
        )

        col4, col5, col6 = st.columns(3)
        ideal_overlap_length = col4.number_input(
            "Preferred overlap length (bp)",
            min_value=1,
            max_value=500,
            value=int(defaults["ideal_overlap_length"]),
            step=1,
            key=f"{prefix}_setting_ideal_length",
        )
        ideal_gc_percent = col5.number_input(
            "Preferred overlap GC (%)",
            min_value=0.0,
            max_value=100.0,
            value=float(defaults["ideal_gc_percent"]),
            step=1.0,
            key=f"{prefix}_setting_ideal_gc",
        )
        max_fragment_size = col6.number_input(
            "Preferred maximum fragment size (bp)",
            min_value=100,
            max_value=1_000_000,
            value=int(defaults["max_fragment_size"]),
            step=100,
            key=f"{prefix}_setting_max_fragment",
            help="Fragments above this size are allowed but receive a score penalty.",
        )

        with st.expander("Advanced design settings", expanded=False):
            st.markdown(
                """
**How overlap-set scoring works**

CombPhage first assigns each candidate overlap a sequence penalty. Lower is better:

`overlap penalty = length penalty + GC penalty`

- For an overlap shorter than the preferred length:  
  `length penalty = missing bp x short-overlap multiplier`
- For an overlap longer than the preferred length:  
  `length penalty = excess bp x long-overlap multiplier`
- `GC penalty = |observed GC fraction - preferred GC fraction| x GC-deviation multiplier`
- When a full filtered overlap is longer than the preferred overlap length, CombPhage tests preferred-length windows and uses the lowest-penalty window for automatic ranking.

The complete overlap-set score is:

`design score = fragment-balance weight x mean physical fragment-size SD`  
`             + overlap-quality weight x mean overlap penalty`  
`             + large-fragment multiplier x total bp above the preferred maximum fragment size`

Fragment-size standard deviation is calculated separately for each parental genome and then averaged across genomes. **Lower design scores rank better.**
                """
            )

            adv1, adv2, adv3 = st.columns(3)
            weight_group_std_dev = adv1.number_input(
                "Fragment-balance weight",
                min_value=0.0,
                value=float(defaults["weight_group_std_dev"]),
                step=0.1,
                key=f"{prefix}_setting_balance_weight",
                help="Weight applied to the mean within-genome fragment-size standard deviation.",
            )
            weight_individual_penalty = adv2.number_input(
                "Overlap-quality weight",
                min_value=0.0,
                value=float(defaults["weight_individual_penalty"]),
                step=1.0,
                key=f"{prefix}_setting_overlap_weight",
                help="Weight applied to the mean sequence penalty of the selected overlaps.",
            )
            large_fragment_multiplier = adv3.number_input(
                "Large-fragment multiplier",
                min_value=0.0,
                value=float(defaults["large_fragment_multiplier"]),
                step=0.1,
                key=f"{prefix}_setting_large_fragment_multiplier",
                help="Penalty per bp exceeding the preferred maximum fragment size.",
            )

            adv4, adv5, adv6 = st.columns(3)
            short_overlap_multiplier = adv4.number_input(
                "Short-overlap multiplier",
                min_value=0.0,
                value=float(defaults["short_overlap_multiplier"]),
                step=1.0,
                key=f"{prefix}_setting_short_multiplier",
                help="Penalty per bp below the preferred overlap length.",
            )
            long_overlap_multiplier = adv5.number_input(
                "Long-overlap multiplier",
                min_value=0.0,
                value=float(defaults["long_overlap_multiplier"]),
                step=0.1,
                key=f"{prefix}_setting_long_multiplier",
                help="Penalty per bp above the preferred overlap length.",
            )
            gc_multiplier = adv6.number_input(
                "GC-deviation multiplier",
                min_value=0.0,
                value=float(defaults["gc_multiplier"]),
                step=10.0,
                key=f"{prefix}_setting_gc_multiplier",
                help=(
                    "Multiplies absolute GC-fraction deviation. For example, "
                    "with 200, a 0.10 (10 percentage-point) deviation adds 20."
                ),
            )

            max_combinations = st.number_input(
                "Maximum overlap combinations to evaluate",
                min_value=1,
                max_value=100_000_000,
                value=int(defaults["max_combinations"]),
                step=10_000,
                key=f"{prefix}_setting_max_combinations",
                help=(
                    "Computational safety limit for exhaustive overlap-set ranking. "
                    "Increasing it can substantially increase runtime."
                ),
            )

    current_settings = {
        "num_fragments": int(num_fragments),
        "topology": topology_label.lower(),
        "top_n": int(top_n),
        "ideal_overlap_length": int(ideal_overlap_length),
        "ideal_gc_percent": float(ideal_gc_percent),
        "max_fragment_size": int(max_fragment_size),
        "min_overlap": int(defaults.get("min_overlap", 1)),
        "weight_group_std_dev": float(weight_group_std_dev),
        "weight_individual_penalty": float(weight_individual_penalty),
        "short_overlap_multiplier": float(short_overlap_multiplier),
        "long_overlap_multiplier": float(long_overlap_multiplier),
        "gc_multiplier": float(gc_multiplier),
        "large_fragment_multiplier": float(large_fragment_multiplier),
        "max_combinations": int(max_combinations),
    }

    st.caption(
        "Re-ranking resets the current working assembly design and any fragment/chimera "
        "outputs derived from it, but does not rerun or modify the upstream CombPhage "
        "overlap analysis."
    )

    if st.button(
        "Rank assembly designs",
        type="primary",
        use_container_width=True,
        key=f"{prefix}_rank_designs",
    ):
        try:
            if outdir.exists():
                shutil.rmtree(outdir)
            initialize_assembly_design(run_dir, current_settings)
            st.success(
                "Assembly designs ranked. Rank 1 is loaded below as the working design."
            )
            st.rerun()
        except Exception as exc:
            st.error(f"Could not rank assembly designs: {exc}")

    context, designs, settings = ensure_assembly_design_state(run_dir)
    if not context or not designs or not settings:
        st.info(
            "Choose the design settings above and click **Rank assembly designs**. "
            "You can repeat this step with different settings without rerunning "
            "the main CombPhage pipeline."
        )
        return

    st.success(
        f"Evaluated {st.session_state[f'{prefix}_total_combinations']:,} "
        f"possible overlap combinations from {len(context['candidates'])} "
        "filtered candidates. Rank 1 is loaded as the editable working design."
    )

    with st.expander("Ranked overlap sets", expanded=True):
        st.dataframe(
            ranked_summary_dataframe(designs, context["genomes"]),
            use_container_width=True,
            hide_index=True,
        )
        st.caption(
            "Lower scores are better. Expand Advanced design settings above to "
            "see the exact scoring terms and coefficients used for this ranking."
        )

    reset_col, info_col = st.columns([1, 3])
    if reset_col.button(
        "Reset to rank 1",
        key=f"{prefix}_reset_rank1",
        use_container_width=True,
    ):
        reset_working_design(
            prefix,
            designs[0],
            context["primary_genome"],
        )
        st.session_state[f"{prefix}_generated_signature"] = None
        st.session_state[f"{prefix}_fragment_analysis_signature"] = None
        st.session_state[f"{prefix}_full_analysis_signature"] = None
        st.rerun()

    info_col.info(
        "Changing a candidate or its working subsequence recalculates the score, "
        "fragment sizes, and overlap coordinates immediately."
    )

    st.subheader("Edit working overlap set")

    candidate_map = context["candidate_map"]
    primary = context["primary_genome"]
    base_ordered = sorted(
        designs[0]["combo"],
        key=lambda overlap: overlap["coords"][primary]["center"],
    )

    selected_ids = []
    selected_candidates = []
    edit_signature = []
    edit_error = False

    def option_label(candidate_id):
        candidate = candidate_map[candidate_id]
        position = candidate["coords"][primary]["center"]
        return (
            f"{candidate_id} | {position:.0f} bp | "
            f"{candidate['length']} bp optimized | "
            f"{candidate['original_length']} bp full | "
            f"GC {candidate['gc'] * 100:.1f}%"
        )

    for index, base_overlap in enumerate(base_ordered):
        reference_position = base_overlap["coords"][primary]["center"]
        options = sorted(
            candidate_map,
            key=lambda candidate_id: abs(
                candidate_map[candidate_id]["coords"][primary]["center"]
                - reference_position
            ),
        )

        widget_key = f"{prefix}_junction_{index}"
        if widget_key not in st.session_state:
            st.session_state[widget_key] = base_overlap["id"]

        selected = st.selectbox(
            f"Junction {index + 1}",
            options=options,
            format_func=option_label,
            key=widget_key,
            help=(
                "Candidates are ordered by distance from this junction's "
                "original rank-1 position."
            ),
        )

        previous_key = f"{prefix}_junction_{index}_previous"
        subseq_key = f"{prefix}_junction_{index}_subsequence"
        selected_candidate = candidate_map[selected]

        if st.session_state.get(previous_key) != selected:
            st.session_state[previous_key] = selected
            # Use the candidate's automatically optimized sequence as the
            # starting point when a different overlap is selected.
            st.session_state[subseq_key] = selected_candidate["seq"]

        if subseq_key not in st.session_state:
            st.session_state[subseq_key] = selected_candidate["seq"]

        selected_ids.append(selected)

        st.markdown("**Full overlap sequence**")
        st.code(selected_candidate["original_seq"], language=None)
        st.caption(
            f"Full filtered overlap: {selected_candidate['original_length']} bp, "
            f"GC {selected_candidate['original_gc'] * 100:.1f}%. "
            "Enter any exact contiguous subsequence below."
        )

        working_subsequence = st.text_area(
            "Working overlap subsequence",
            key=subseq_key,
            height=90,
            help=(
                "The sequence must occur exactly and contiguously within the "
                "full overlap shown above. Coordinates for every genome are "
                "updated automatically."
            ),
        )

        try:
            edited_candidate = candidate_from_subsequence(
                selected_candidate,
                context["genome_lengths"],
                subsequence=working_subsequence,
                ideal_overlap_length=context["ideal_overlap_length"],
                ideal_gc=context["ideal_gc"],
                min_length=int(settings.get("min_overlap", 1)),
                short_overlap_multiplier=float(settings["short_overlap_multiplier"]),
                long_overlap_multiplier=float(settings["long_overlap_multiplier"]),
                gc_multiplier=float(settings["gc_multiplier"]),
            )
        except Exception as exc:
            st.error(str(exc))
            edit_error = True
            continue

        st.caption(
            f"Selected subsequence: {edited_candidate['length']} bp, "
            f"GC {edited_candidate['gc'] * 100:.1f}%, "
            f"offset {edited_candidate['left_trim']} bp from the left and "
            f"{edited_candidate['right_trim']} bp from the right of the full overlap."
        )

        selected_candidates.append(edited_candidate)
        edit_signature.append(
            (
                selected,
                edited_candidate["seq"],
                int(edited_candidate["left_trim"]),
                int(edited_candidate["right_trim"]),
            )
        )

    has_duplicates = len(selected_ids) != len(set(selected_ids))

    if has_duplicates:
        st.error(
            "The same overlap has been selected for more than one junction. "
            "Choose a unique overlap for each junction."
        )
        current_design = None
    elif edit_error or len(selected_candidates) != len(base_ordered):
        current_design = None
    else:
        try:
            current_design = design_from_candidates(
                context,
                selected_candidates,
                topology=settings["topology"],
                max_fragment_size=settings["max_fragment_size"],
                weight_group_std_dev=float(settings["weight_group_std_dev"]),
                weight_individual_penalty=float(settings["weight_individual_penalty"]),
                large_fragment_multiplier=float(settings["large_fragment_multiplier"]),
            )
        except Exception as exc:
            st.error(f"Could not evaluate the edited design: {exc}")
            current_design = None

    if current_design is None:
        return

    score_cols = st.columns(4)
    score_cols[0].metric(
        "Current score",
        f"{current_design['score']:.2f}",
    )
    score_cols[1].metric(
        "Rank-1 score",
        f"{designs[0]['score']:.2f}",
    )
    score_cols[2].metric(
        "Mean fragment-size SD",
        f"{current_design['group_std_dev']:.1f} bp",
    )
    score_cols[3].metric(
        "Overlap penalty",
        f"{current_design['average_overlap_penalty']:.2f}",
    )

    junction_table, fragment_table = current_design_tables(
        current_design,
        context,
    )

    display_left, display_right = st.columns([1.25, 1])

    with display_left:
        st.markdown("**Working junctions in genomic order**")
        st.dataframe(
            junction_table,
            use_container_width=True,
            hide_index=True,
        )

    with display_right:
        st.markdown("**Predicted physical fragment sizes (bp; overlaps included)**")
        st.caption(
            "These are the exact fragment lengths used for FASTA extraction. "
            "Each shared assembly overlap is retained on both adjacent fragments."
        )
        st.dataframe(
            fragment_table,
            use_container_width=True,
            hide_index=True,
        )

    current_signature = tuple(edit_signature)

    fragment_signature = st.session_state.get(
        f"{prefix}_fragment_analysis_signature"
    )
    full_signature = st.session_state.get(
        f"{prefix}_full_analysis_signature"
    )

    if (
        fragment_signature is not None
        and current_signature != tuple(fragment_signature)
    ):
        st.warning(
            "The working design has changed since the fragment analysis was "
            "generated. Re-analyze the fragments to evaluate the current design."
        )

    if (
        full_signature is not None
        and current_signature != tuple(full_signature)
    ):
        st.warning(
            "The working design has changed since the full chimera library was "
            "generated. The existing library/analysis represents an older design."
        )

    st.subheader("1. Analyze current fragment design")
    st.write(
        "Generate only the parental versions of each fragment and analyze "
        "fragment nucleotide diversity, homologous CDS diversity, and the most "
        "variable CDSs. This is intended to be rerun while refining junctions."
    )

    analysis_options = st.columns(2)
    top_variable_cds = analysis_options[0].number_input(
        "Variable CDSs to show per fragment",
        min_value=1,
        max_value=25,
        value=DEFAULT_TOP_VARIABLE_CDS,
        step=1,
        key=f"{prefix}_top_variable_cds",
    )
    save_alignments = analysis_options[1].checkbox(
        "Save pairwise fragment alignments",
        value=False,
        key=f"{prefix}_save_fragment_alignments",
    )

    if st.button(
        "Generate & analyze current fragments",
        type="primary",
        use_container_width=True,
        key=f"{prefix}_analyze_fragments",
    ):
        try:
            working_dir = outdir / "working_design"
            materialize_design(
                current_design,
                context["genome_sequences"],
                context["genomes"],
                settings["topology"],
                working_dir,
                generate_chimeras=False,
                clean_outdir=True,
            )
            write_single_design(
                current_design,
                context["genomes"],
                working_dir / "working_overlaps.tsv",
            )

            analyze_fragment_design(
                materialized_dir=working_dir,
                prepped_dir=run_dir / "prepped_genbanks",
                cluster_membership=(
                    run_dir / "feature_overlaps" / "cluster_membership.tsv"
                ),
                outdir=working_dir / "analysis",
                save_alignments=bool(save_alignments),
                top_variable_cds=int(top_variable_cds),
            )

            st.session_state[
                f"{prefix}_fragment_analysis_signature"
            ] = current_signature
            st.success(
                "Current fragments generated and analyzed. Edit the junctions "
                "and run this step again if you want to compare another design."
            )
        except Exception as exc:
            st.error(f"Could not analyze current fragment design: {exc}")

    show_fragment_analysis_outputs(run_dir, prefix)

    st.subheader("2. Generate & analyze full assembly library")
    st.write(
        "Once the fragment design is satisfactory, generate every non-parental "
        "fragment combination and calculate assembly-vs-parent and, when the "
        "library is small enough, all-vs-all assembly identities."
    )

    generate_col1, generate_col2, generate_col3 = st.columns(3)
    chimera_prefix = generate_col1.text_input(
        "Chimera prefix",
        value="Chimera",
        key=f"{prefix}_chimera_prefix",
    )
    max_chimeras = generate_col2.number_input(
        "Maximum chimeras to generate",
        min_value=1,
        max_value=1_000_000,
        value=DEFAULT_MAX_CHIMERAS,
        step=100,
        key=f"{prefix}_max_chimeras",
    )
    max_pairwise_assemblies = generate_col3.number_input(
        "Max assemblies for all-vs-all analysis",
        min_value=2,
        max_value=10_000,
        value=DEFAULT_MAX_PAIRWISE_ASSEMBLIES,
        step=25,
        key=f"{prefix}_max_pairwise_assemblies",
    )

    possible_chimeras = (
        len(context["genomes"]) ** settings["num_fragments"]
        - len(context["genomes"])
    )
    st.caption(
        f"This design produces {possible_chimeras:,} non-parental fragment "
        "combinations."
    )

    full_disabled = possible_chimeras > int(max_chimeras)
    if full_disabled:
        st.warning(
            f"This exceeds the current limit of {int(max_chimeras):,} chimeras. "
            "Increase the limit to enable full-library generation."
        )

    if st.button(
        "Generate & analyze full assembly library",
        type="primary",
        use_container_width=True,
        key=f"{prefix}_generate_and_analyze_full",
        disabled=full_disabled,
    ):
        try:
            selected_dir = outdir / "selected_design"
            materialize_design(
                current_design,
                context["genome_sequences"],
                context["genomes"],
                settings["topology"],
                selected_dir,
                chimera_prefix=chimera_prefix,
                generate_chimeras=True,
                max_chimeras=int(max_chimeras),
                clean_outdir=True,
            )
            write_single_design(
                current_design,
                context["genomes"],
                selected_dir / "selected_overlaps.tsv",
            )

            analyze_materialized_design(
                materialized_dir=selected_dir,
                prepped_dir=run_dir / "prepped_genbanks",
                cluster_membership=(
                    run_dir / "feature_overlaps" / "cluster_membership.tsv"
                ),
                outdir=selected_dir / "analysis",
                save_alignments=bool(save_alignments),
                top_variable_cds=int(top_variable_cds),
                max_pairwise_assemblies=int(max_pairwise_assemblies),
            )

            st.session_state[f"{prefix}_generated_signature"] = current_signature
            st.session_state[f"{prefix}_full_analysis_signature"] = current_signature
            st.success("Full chimera library generated and analyzed.")
        except Exception as exc:
            st.error(f"Could not generate/analyze the full assembly library: {exc}")

    show_full_design_analysis_outputs(run_dir, prefix)


def show_complete_run_download(run_dir: Path):
    archive_path = archive_directory(
        run_dir,
        f"{run_dir.name}_results",
    )

    st.download_button(
        "Download complete CombPhage run (.zip)",
        data=archive_path.read_bytes(),
        file_name=archive_path.name,
        mime="application/zip",
        key=f"complete_{run_dir.name}",
    )


st.set_page_config(page_title="CombPhage", layout="wide")
st.title("CombPhage")
st.caption(
    "Identify conserved assembly overlaps and design interchangeable "
    "bacteriophage genome fragments."
)

st.subheader("Genome input")
input_mode = st.radio(
    "Input method",
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

    email = st.text_input(
        "Email for NCBI Entrez",
        placeholder="email@example.com",
    )

st.divider()

with st.expander("Overlap analysis settings", expanded=True):
    overlap_col1, overlap_col2, overlap_col3 = st.columns(3)

    min_overlap = overlap_col1.number_input(
        "Minimum overlap size (bp)",
        min_value=1,
        max_value=500,
        value=25,
        step=1,
    )

    max_homopolymer = overlap_col2.number_input(
        "Maximum allowed homopolymer length",
        min_value=1,
        max_value=20,
        value=4,
        step=1,
    )

    generate_viz = overlap_col3.checkbox(
        "Generate LoVis4u visualization",
        value=True,
    )

with st.expander("Run settings", expanded=False):
    runs_dir_text = st.text_input(
        "Runs directory",
        value=str(DEFAULT_RUNS_DIR),
    )


run_name = st.text_input(
    "Run name",
    value="combphage_" + datetime.now().strftime("%Y%m%d_%H%M%S"),
)


if st.button(
    "Run CombPhage",
    type="primary",
    use_container_width=True,
):
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
            st.error(
                f"Run directory already exists: {run_dir}\n"
                "Choose a different run name."
            )
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
                        f"{icons[stage_states[stage]]} "
                        f"{STAGE_LABELS[stage]}"
                        for stage in STAGES
                    )
                )

            def callback(stage, state, message):
                stage_states[stage] = state
                status.write(message)
                redraw()

                done = sum(
                    state_name in {"complete", "skipped"}
                    for state_name in stage_states.values()
                )
                progress.progress(done / len(STAGES))

            redraw()

            try:
                pipeline = CombPhagePipeline(
                    REPO_ROOT,
                    run_dir,
                    callback,
                )

                if input_mode == "Upload GenBank files":
                    incoming = run_dir / "_uploaded_genbanks"
                    save_uploads(uploaded_files, incoming)
                    pipeline.use_local_genbanks(incoming)
                else:
                    pipeline.fetch_genbanks(
                        parse_accessions(accession_text),
                        email,
                    )

                pipeline.run_all(
                    min_size=int(min_overlap),
                    max_homo_size=int(max_homopolymer),
                    visualize=generate_viz,
                )

                progress.progress(1.0)
                status.success(
                    "CombPhage overlap analysis complete. Configure assembly design below when ready."
                )
                st.session_state["last_combphage_run"] = str(run_dir)

            except CombPhagePipelineError as exc:
                st.error(str(exc))

                if exc.log_file and exc.log_file.exists():
                    with st.expander(
                        f"{exc.stage} log",
                        expanded=True,
                    ):
                        st.code(
                            exc.log_file.read_text(
                                encoding="utf-8",
                                errors="replace",
                            )[-12000:]
                        )

            except Exception as exc:
                if run_dir.exists():
                    st.session_state["last_combphage_run"] = str(run_dir)
                st.exception(exc)


last_run = st.session_state.get("last_combphage_run")

if last_run:
    run_path = Path(last_run)

    if run_path.exists():
        st.divider()
        show_results(run_path)
        show_assembly_design(run_path)
        st.divider()
        show_complete_run_download(run_path)
