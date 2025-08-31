# MethanoPipeline: Installation Guide

This guide helps you set up **MethanoPipeline** on a fresh machine. Follow these steps to get the pipeline running with minimal hassle.

---

## Prerequisites

Install these before cloning the repository:

* **R** (≥ 4.2.0): [cran.r-project.org](https://cran.r-project.org/)
* **Conda** (Miniconda recommended): [docs.conda.io/miniconda](https://docs.conda.io/en/latest/miniconda.html)
* **Git**: [git-scm.com/downloads](https://git-scm.com/downloads)

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/timotheqds/Methano-Project-16s-Data-Analysis.git
cd Methano-Project-16s-Data-Analysis

# 2. Create and activate conda environment
conda env create -f environment.yml
conda activate methanopipeline

# 3. Install required R packages
Rscript scripts/install_r_packages.R

# 4. (Optional) Download reference databases for taxonomy assignment
GTDB_r226-modApril2025-RData and silva_species_assignment_v138-1.fa.gz

# 5. Verify installation
Rscript scripts/verify_installation.R
```

---

## Data Setup

1. Place FASTQ files in `data/raw/`

   * Naming: `SAMPLEID_R1.fastq.gz`, `SAMPLEID_R2.fastq.gz`
2. Create your metadata file using `data/metadata_template.csv` as a guide.

---

## Running the Pipeline

Run run an individual step (useful for testing):

```bash
Rscript scripts/steps/step1_setup_validation.R
```

---

## Output Structure

```
data/
├── raw/          # Input FASTQ files
├── processed/    # Intermediate outputs
├── results/      # Final outputs
└── step_outputs/ # Step-by-step results

logs/
└── steps/        # Log files
```

---

## Customization

Edit the following for your dataset:

* `data/metadata_template.csv` → sample metadata
* `scripts/steps/step1_setup_validation.R` → expected sample IDs & primer sequences

See [README.md](../README.md) for full details.

---

## Troubleshooting

If something breaks:

* Check logs in `logs/steps/`
* Refer to `docs/troubleshooting.md` for common issues
* Ensure dependencies and database paths are correct

---
