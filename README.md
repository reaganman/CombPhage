# CombPhage

CombPhage is a comparative bacteriophage genome analysis pipeline designed to identify conserved nucleotide overlaps that can be used as assembly junctions for constructing interchangeable or recombinable phage genomes.

The pipeline compares related phage genomes, groups homologous coding sequences, identifies conserved nucleotide sequence surrounding shared CDS boundaries, filters those candidate overlaps for sequence properties that may complicate DNA assembly, and generates genome-level visualizations highlighting candidate assembly junctions.

CombPhage is intended to support workflows where related phage genomes are divided into interchangeable genomic fragments using conserved sequence overlaps.

---

## What CombPhage Does

CombPhage takes a set of related annotated phage genomes and performs the following workflow:

1. Retrieve GenBank records from NCBI or accept user-provided GenBank files.
2. Normalize genome and CDS identifiers across all input genomes.
3. Extract nucleotide sequences, protein sequences, and CDS metadata.
4. Cluster homologous proteins using MMseqs2.
5. Identify clusters containing corresponding CDSs across the input genomes.
6. Resolve ambiguous clusters containing multiple related CDSs within a genome.
7. Search nucleotide sequence surrounding homologous CDS boundaries for conserved overlaps.
8. Filter candidate overlaps based on sequence properties and minimum overlap length.
9. Generate tables describing candidate assembly junctions.
10. Visualize the genomes and highlight CDSs associated with usable overlaps using LoVis4u.

A simplified workflow is:

```text
Annotated phage genomes
        |
        v
Normalize annotations
        |
        v
Protein clustering
        |
        v
Resolve homologous CDS groups
        |
        v
Find conserved CDS-boundary sequence
        |
        v
Filter candidate overlaps
        |
        v
Candidate interchangeable
genome assembly junctions
        |
        v
Genome visualization
```

---

# Installation

CombPhage is designed to run in a Conda environment.

Clone the repository:

```bash
git clone https://github.com/reaganman/CombPhage.git
cd CombPhage
```

Create the environment:

```bash
conda env create -f environment.yml
```

Activate it:

```bash
conda activate combphage
```

CombPhage currently relies on several external tools and Python packages, including:

- Python
- Biopython
- pandas
- MMseqs2
- LoVis4u
- Streamlit

LoVis4u may require its Linux setup step after installation:

```bash
lovis4u --linux
```

This should normally only need to be performed once for the environment.

---

# Running CombPhage

CombPhage can currently be run either through the command line or through the Streamlit graphical interface.

---

## Command-Line Workflow

The main workflow consists of five steps.

### 1. Obtain GenBank files

If the genomes are available through NCBI, they can be downloaded using accession numbers:

```bash
python Scripts/fetch_genbanks.py \
    -a MZ501081.1 MZ501078.1 V01146.1 \
    -o genbanks \
    --email your@email.com
```

The email address is required by NCBI Entrez.

Alternatively, existing annotated GenBank files can be placed directly in a directory and used without this step.

For example:

```text
genbanks/
├── MZ501081.1.gbk
├── MZ501078.1.gbk
└── V01146.1.gbk
```

---

### 2. Prepare and normalize the genomes

Run:

```bash
python Scripts/prepare_genbanks.py \
    -i genbanks/ \
    -o prepped_genbanks
```

This step standardizes CDS identifiers and creates CombPhage's internal input files.

Typical outputs include:

```text
prepped_genbanks/
├── genomes.fasta
├── proteins.faa
├── cds.tsv
└── genbanks/
```

The normalized CDS identifiers follow the general format:

```text
GENOME_CDS_####
```

For example:

```text
MZ501078.1_CDS_0042
```

The same identifier is used consistently across the prepared FASTA files, CDS table, normalized GenBank files, clustering results, and downstream analyses.

---

### 3. Identify conserved CDS-boundary overlaps

Run:

```bash
python Scripts/find_overlaps.py \
    --input-dir prepped_genbanks/ \
    --outdir feature_overlaps
```

This step:

- clusters related proteins using MMseqs2,
- identifies CDS families shared across genomes,
- distinguishes direct and ambiguous clusters,
- attempts to resolve ambiguous clusters using pairwise protein alignment coverage and identity,
- and searches nucleotide sequence around shared CDS starts and stops for conserved sequence.

The primary candidate-overlap output is:

```text
feature_overlaps/cluster_overlap_summary.tsv
```

Additional diagnostic files describe cluster membership and cases where homologous CDS assignments were incomplete or ambiguous.

---

### 4. Filter candidate overlaps

Run:

```bash
python Scripts/filter_overlaps.py \
    --overlaps feature_overlaps/cluster_overlap_summary.tsv \
    --min_size 25
```

The `--min_size` argument specifies the minimum acceptable overlap length in base pairs.

For example:

```bash
--min_size 25
```

requires candidate overlaps to contain at least 25 bp after filtering.

The main output is:

```text
feature_overlaps/filtered_overlaps.tsv
```

A detailed filtering report is also generated:

```text
feature_overlaps/overlap_filter_report.tsv
```

The filtering step evaluates candidate overlaps for properties that may complicate DNA assembly and can trim unsuitable sequence from otherwise usable overlaps.

---

### 5. Visualize candidate overlaps

Run:

```bash
python Scripts/visualize.py \
    --genbanks prepped_genbanks/genbanks/ \
    --overlaps-dir feature_overlaps
```

This generates a comparative genome visualization using LoVis4u.

CDSs are color-coded based on the location of suitable assembly overlaps:

- **Green** — suitable overlap at the CDS start
- **Red** — suitable overlap at the CDS stop
- **Purple** — suitable overlaps at both the CDS start and stop
- **Grey** — no suitable filtered overlap

Only CDSs associated with suitable overlaps are labeled, and their labels correspond to the CombPhage cluster identifier.

The resulting visualization can be used to quickly identify candidate boundaries for dividing related phage genomes into interchangeable fragments.

---

# Running the GUI

CombPhage also includes a Streamlit interface.

From the repository root:

```bash
conda activate combphage
streamlit run app.py
```

The interface allows genomes to be provided either by:

- uploading annotated GenBank files, or
- entering NCBI nucleotide accession numbers.

The GUI then runs the same underlying CombPhage pipeline used by the command-line workflow.

Analysis settings currently include parameters such as:

- minimum overlap length,
- maximum allowed homopolymer length,
- and whether to generate the LoVis4u visualization.

Each analysis is stored in a separate run directory.

---

# Output Files

The main analysis directory is typically:

```text
feature_overlaps/
```

Important files include the following.

## `cluster_membership.tsv`

Lists the CDSs assigned to each final homologous CDS cluster.

This is useful for determining which genes from different genomes are being treated as corresponding features.

Typical information includes:

```text
cluster
gene
genome
resolution
source_mmseqs_representative
```

Clusters may originate directly from MMseqs2 or may have been resolved from a broader ambiguous protein family.

---

## `cluster_overlap_summary.tsv`

Contains the conserved nucleotide overlaps identified around CDS starts and stops before sequence-quality filtering.

Each row represents a candidate overlap associated with a homologous CDS cluster.

Important fields include:

```text
cluster
region
reference_genome
reference_gene
ov_seq
overlap_length
```

Genome-specific overlap coordinates are also included.

The `region` column identifies whether the conserved sequence is associated with the:

```text
start
```

or:

```text
stop
```

of the corresponding CDS.

---

## `filtered_overlaps.tsv`

This is the primary result file for selecting candidate assembly junctions.

It contains overlaps that remain suitable after filtering and trimming.

Important information includes:

- cluster identifier,
- CDS boundary,
- overlap sequence,
- overlap length,
- reference CDS,
- genome-specific coordinates,
- and cluster-resolution metadata.

For most downstream experimental design, this is the most useful table to begin with.

---

## `overlap_filter_report.tsv`

Contains a record of every candidate overlap evaluated by the filtering step, including rejected overlaps.

The report can be used to determine why a candidate was:

- accepted,
- trimmed,
- or rejected.

It also preserves the original sequence and original overlap length so that filtering decisions can be examined.

Typical fields include:

```text
cluster
region
status
reason
filter_actions
original_ov_seq
final_ov_seq
original_overlap_length
final_overlap_length
left_trim
right_trim
```

This file is particularly useful when an expected overlap does not appear in `filtered_overlaps.tsv`.

---

## `ambiguous_clusters.tsv`

Contains broad MMseqs2 protein clusters where more than one CDS from the same genome was assigned to the same protein family.

These may represent:

- paralogous genes,
- overlapping alternative annotations,
- split or fused CDS annotations,
- or closely related proteins within the same genome.

CombPhage attempts to resolve these families before using them for overlap identification.

---

## `resolved_clusters.tsv`

Describes CDS assignments selected when an ambiguous MMseqs2 cluster can be separated into one or more consistent cross-genome groups.

Resolution is based on pairwise protein similarity, including balanced sequence coverage and identity.

---

## `unresolved_ambiguous_members.tsv`

Contains CDSs from ambiguous protein families that could not be confidently assigned to a resolved cross-genome cluster.

These are excluded from downstream overlap identification.

---

## `incomplete_clusters.tsv`

Contains protein clusters that do not contain CDSs from enough of the input genomes to be used as shared interchangeable boundaries.

These clusters can still be biologically informative but are not used as fully conserved assembly junctions.

---

## `unclustered_cds.tsv`

Contains CDSs that were not assigned to an MMseqs2 cluster.

---

# Interpreting the Results

CombPhage is designed to help identify genomic boundaries that could serve as interchangeable assembly junctions.

A useful overlap should generally meet two requirements:

1. The corresponding CDS boundary is homologous across the genomes being compared.
2. The nucleotide sequence surrounding that boundary is sufficiently conserved and suitable for DNA assembly.

The most direct starting point is:

```text
filtered_overlaps.tsv
```

Each retained row represents a candidate conserved sequence that could potentially be used as an overlap between neighboring genome fragments.

For example, if a cluster has:

```text
cluster12    start
```

then a conserved sequence exists around the start of the CDSs belonging to `cluster12`.

If the same cluster also has:

```text
cluster12    stop
```

then suitable conserved overlaps occur on both sides of that homologous CDS.

These junctions can be used to define fragment boundaries.

---

# Designing Interchangeable Genome Fragments

The genome visualization is useful for identifying combinations of overlaps that partition related genomes into corresponding fragments.

For example:

```text
Genome A
------[cluster3]-----------[cluster8]-----------[cluster14]------

Genome B
------[cluster3]-----------[cluster8]-----------[cluster14]------

Genome C
------[cluster3]-----------[cluster8]-----------[cluster14]------
```

If suitable overlaps occur at each of these conserved positions, the genomes could conceptually be divided into:

```text
Fragment 1
cluster3 → cluster8

Fragment 2
cluster8 → cluster14

Fragment 3
cluster14 → next junction
```

Corresponding fragments from related genomes can then be treated as candidate interchangeable units.

CombPhage identifies conserved candidate junctions; experimental feasibility should still be evaluated based on the intended DNA assembly method, genome biology, and experimental design.

---

# Choosing Between Multiple Candidate Overlaps

When several overlaps are available near a desired genome boundary, useful factors to consider include:

- overlap length,
- sequence composition,
- whether sequence filtering required trimming,
- whether the cluster was directly identified or computationally resolved,
- whether the same junction is present in all genomes of interest,
- and the biological location of the boundary.

The `overlap_filter_report.tsv` file can help evaluate why particular overlaps were retained or rejected.

The visualization can then be used to determine whether the retained overlaps produce a useful partitioning of the genomes.

---

# Example Workflow

A complete command-line analysis might look like:

```bash
conda activate combphage

python Scripts/fetch_genbanks.py \
    -a MZ501081.1 MZ501078.1 V01146.1 \
    -o genbanks \
    --email your@email.com

python Scripts/prepare_genbanks.py \
    -i genbanks/ \
    -o prepped_genbanks

python Scripts/find_overlaps.py \
    --input-dir prepped_genbanks/ \
    --outdir feature_overlaps

python Scripts/filter_overlaps.py \
    --overlaps feature_overlaps/cluster_overlap_summary.tsv \
    --min_size 25

python Scripts/visualize.py \
    --genbanks prepped_genbanks/genbanks/ \
    --overlaps-dir feature_overlaps
```

The primary outputs to inspect afterward are:

```text
feature_overlaps/filtered_overlaps.tsv
feature_overlaps/overlap_filter_report.tsv
feature_overlaps/cluster_membership.tsv
feature_overlaps/lovis4u_overlaps/
```

---

# Current Development Status

CombPhage is under active development.

The current implementation focuses on:

- related bacteriophage genomes,
- annotated GenBank input,
- homologous CDS-based genome comparison,
- conserved CDS-boundary overlap discovery,
- protein-family ambiguity resolution,
- overlap sequence filtering,
- and comparative genome visualization.

Future development may include expanded interactive visualization and additional support for downstream interchangeable genome design.

---

# Citation

A formal citation for CombPhage will be added as the software and associated methodology are developed.

---

# License

CombPhage is released under the MIT License. See the LICENSE file for details.
