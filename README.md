# MethanoPipeline: 16S rRNA Analysis of Tree Soil, Bark, and Wood

A customised bioinformatics pipeline for analysing **16S rRNA amplicon sequencing data** from tree-associated environments (soil, bark, and wood). Built upon **SimpleMetaPipeline**, optimized for **environmental microbiome data**, and designed for reproducibility, quality control, and downstream analysis.

---

## Prerequisites

### Software

* **Conda** (Miniconda/Anaconda)
* **R (v4.2+)** with packages: `dada2`, `phyloseq`, `DECIPHER`, `lulu`
* **Git**

### Tools (via Conda)

* **cutadapt** – primer/adaptor removal
* **fastqc** / **multiqc** – quality control
* **swarm** – OTU clustering
* **blast+** – fallback taxonomy assignment and LULU matchlist support

### Databases

* **GTDB** (primary taxonomy)
* **SILVA** (fallback taxonomy)
* **BLAST** (final rescue for unclassified sequences)

---

## Data Requirements & Structure

* **Paired-end FASTQ files**: `SAMPLE_R1.fastq.gz`, `SAMPLE_R2.fastq.gz`
* **Metadata CSV**: `data/metadata.csv` with columns
  `sample-id, environment, tree-species, collection-date, location`

**Example:**

```csv
sample-id,environment,tree-species,collection-date,location
A6S,soil,oak,2024-07-01,Plot1
B6B,bark,pine,2024-07-01,Plot2
B10W,wood,beech,2024-07-02,Plot3
```

```
project/
├─ data/
│  ├─ raw/                  # FASTQ files
│  ├─ processed/            # Intermediate files
│  ├─ results/              # Final outputs
│  └─ step_outputs/         # Per-step results
└─ logs/
   └─ steps/                # Execution logs
```

---

## Pipeline Scripts (at a glance)

| Step | Script                  | Purpose                                                      |
| ---- | ----------------------- | ------------------------------------------------------------ |
| 1    | `0.1_setup.R`           | Environment setup, metadata validation, paths                |
| 2    | `0.2_quality_control.R` | Primer removal (cutadapt), QC reports (FastQC/MultiQC)       |
| 3    | `0.3_filtering.R`       | Quality filtering, trimming, dereplication (DADA2)           |
| 4    | `0.4_sequence_table.R`  | Error learning → ASV inference, merge pairs, chimera removal |
| 5    | `0.5_lulu_curation.R`   | LULU curation (stringent + relaxed)                          |
| 6    | `0.6_otu_clustering.R`  | Swarm OTU clustering (optional, alongside ASVs)              |
| 7    | `0.7_taxonomy.R`        | Taxonomy with IDTAXA (GTDB) + SILVA fallback + BLAST rescue  |

---

## Quick Start

```bash
# Clone repository
git clone https://github.com/timotheqds/Methano-Project-16s-Data-Analysis.git
cd Methano-Project-16s-Data-Analysis

# Set up environment
conda env create -f environment.yml
conda activate methanopipeline
Rscript scripts/install_r_packages.R

# Configure (sample names, primers) in scripts/0.1_setup.R, then run:
Rscript scripts/0.1_setup.R
Rscript scripts/0.2_quality_control.R
Rscript scripts/0.3_filtering.R
Rscript scripts/0.4_sequence_table.R
Rscript scripts/0.5_lulu_curation.R
Rscript scripts/0.6_otu_clustering.R
Rscript scripts/0.7_taxonomy.R
```

**Example configuration (`scripts/0.1_setup.R`)**

```r
expected_samples <- c("A6S", "B6B", "B10W")  # example
FWD <- "AACMGGATTAGATACCCKG"   # 799F
REV <- "ACGTCRTCCMCACCTTCCTC" # 1193R
```

---

## Step-by-step Explanations

### Step 1 — Environment Setup & Data Preparation (`0.1_setup.R`)

**Goal:** Make sure everything is in the right place and consistent before heavy processing.

* **Inputs:** `data/raw/*.fastq.gz`, `data/metadata.csv`, primer sequences.

* **Actions:**

  * Validate FASTQ naming (e.g., `_R1/_R2` pairs) and file integrity.
  * Parse metadata; verify all `sample-id`s are present in FASTQs and vice versa.
  * Create/clean output directories under `data/processed`, `data/results`, `data/step_outputs`, and `logs/steps`.
  * Initialize logging and write session info (R/Conda versions) for reproducibility.

* **Outputs:** `logs/steps/00_setup.log`, `data/step_outputs/setup/sample_map.tsv`.

* **Tips:** Keep `expected_samples` explicit to catch missing or extra files early.

---

### Step 2 — Quality Control & Primer Removal (`0.2_quality_control.R`)

**Goal:** Remove primers/adapters and assess raw read quality.

* **Inputs:** Paired FASTQs from `data/raw/`.
* **Actions:**

  * **Cutadapt**: remove forward/reverse primers; optionally anchor to read starts; allow small mismatch rate (e.g., 0.1).
  * **FastQC**: generate per-sample QC metrics (per-base quality, GC, length dist.).
  * **MultiQC**: aggregate reports across samples.
* **Outputs:**

  * Primer-trimmed FASTQs in `data/processed/trimmed/`.
  * QC reports in `data/results/quality_control/` (FastQC HTMLs, MultiQC summary).
* **Watch:** Over-trimming (short reads), residual primers (non-anchored matches), and adapter content flags.

---

### Step 3 — Read Filtering & Trimming (`0.3_filtering.R`)

**Goal:** Improve data quality before denoising.

* **Actions (DADA2 `filterAndTrim`)**:

  * Trim low-quality tails by examining per-base quality plots.
  * Set reasonable thresholds (examples): `maxEE` \~ 2 (R1) / 4 (R2); `truncQ` \~ 2–11; remove `Ns`; `rm.phix=TRUE`.
  * Dereplicate identical sequences to speed inference.
* **Optional abundance screens** (applied later if preferred):

  * ≥ 5 reads per sample **and** ≥ 10 total reads across all samples **and** ≥ 0.1% relative abundance.
* **Outputs:** Filtered FASTQs in `data/processed/filtered/`, summary `data/step_outputs/filtering/read_retention.tsv`.
* **Watch:** Don’t truncate so aggressively that read-pair overlap becomes impossible (see next step).

---

### Step 4 — Sequence Table Construction & Denoising (`0.4_sequence_table.R`)

**Goal:** Infer high-resolution ASVs and remove artifacts.

* **Actions:**

  * Learn error models per direction (`learnErrors`).
  * Infer ASVs per sample (`dada`).
  * **Merge paired reads** where possible (`mergePairs`).
  * Remove chimeras (`removeBimeraDenovo`).
  * Construct the **ASV table** (samples × ASVs).
* **Automatic fallback:** If the **merge success rate < 10%**, the pipeline switches to **forward-reads-only**. This typically retains \~97% of reads and captures V5–V6 (\~250 bp) with adequate taxonomic resolution.
* **Outputs:**

  * `data/processed/asv/seqtab.rds`, `data/results/ASV_table.csv` (or `ASV_table_curated.csv` after LULU), chimera/retention stats.
* **Watch:** Extremely low overlap with 2×250 bp runs of \~400+ bp amplicons; consider 2×300 bp, alternative primer sets (e.g., 515F–806R), or long-read sequencing.

---

### Step 5 — LULU Post-clustering Curation (`0.5_lulu_curation.R`)

**Goal:** Reduce spurious ASVs by exploiting co-occurrence and sequence similarity.

* **Actions:**

  * Build a **matchlist** (via BLAST or vsearch) of close ASV sequences.
  * Run **LULU** first with conservative thresholds, then with relaxed thresholds.
  * Compare curated vs. uncurated tables; track collapsed ASVs.
* **Outputs:** `data/results/ASV_table_curated.csv`, curation report, matchlist artifact files.
* **Watch:** Over-aggressive curation can merge genuine ecological variants—inspect per-sample effects.

---

### Step 6 — OTU Clustering with Swarm (`0.6_otu_clustering.R`)

**Goal:** Complement ASV view with OTU-level groupings.

* **Actions:**

  * Cluster ASVs using **Swarm** (e.g., d=1) to produce OTUs.
  * Generate **OTU table** and mapping from ASV→OTU.
  * Compare alpha/beta diversity from ASVs vs. OTUs.
* **Outputs:** `data/results/OTU_table.csv`, `data/step_outputs/swarm/`.
* **Watch:** Swarm parameters (d) strongly affect richness; document chosen value in logs.

---

### Step 7 — Taxonomic Assignment (`0.7_taxonomy.R`)

**Goal:** Assign taxonomy to ASVs/OTUs using robust references.

* **Actions:**

  * **IDTAXA (DECIPHER)** with **GTDB** training set as primary.
  * **SILVA fallback** for sequences unclassified with GTDB.
  * **BLAST+ rescue** for any sequences still unclassified.
  * Report confidence/support per rank (domain→species), and success rates.
* **Outputs:** `data/results/taxonomy_assignment_gtdb.csv` (plus SILVA- and BLAST-enhanced table), per-rank summaries.
* **Watch:** Keep GTDB/SILVA versions pinned; note bootstrap thresholds used.

---

## Required Configuration

Edit `scripts/0.1_setup.R`:

```r
expected_samples <- c("SAMPLE1", "SAMPLE2", "SAMPLE3")
FWD <- "AACMGGATTAGATACCCKG"  # 799F primer
REV <- "ACGTCRTCCMCACCTTCCTC"  # 1193R primer
```

Optional abundance filters (apply during/after Step 4):

* **Per-sample minimum:** ≥ 5 reads
* **Global minimum:** ≥ 10 reads total
* **Relative abundance:** ≥ 0.1%

---

## Key Outputs (summary)

* `data/results/ASV_table_curated.csv` – curated ASV abundance table
* `data/results/OTU_table.csv` – OTU abundance table (if clustering run)
* `data/results/taxonomy_assignment_gtdb.csv` – taxonomy results (with SILVA + BLAST fallback)
* `logs/steps/` – detailed logs for reproducibility

---

## Notes & Best Practices

* **Automatic forward-only fallback** if paired-end merging success < 10%, typically retaining \~97% of reads (V5–V6, \~250 bp).
* **Merging challenges:** For \~400 bp amplicons with 2×250 bp reads, consider **2×300 bp**, alternative **515F–806R** primer set, or **long-read** platforms.
* **Reproducibility:** Pin database/tool versions; commit `sessionInfo()` and `conda env export` snapshots.
* **Resumability:** Each step writes checkpoints; you can rerun from any step without re-processing earlier outputs.

---

## Documentation

* `INSTALL.md` – Installation instructions
* `CONFIGURATION.md` – Customization & setup
* `docs/troubleshooting.md` – Common issues and fixes

---

## Attribution

Built on concepts from **SimpleMetaPipeline** with customizations for tree-associated microbiomes (soil, bark, wood).
