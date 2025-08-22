# MethanoPipeline: 16S rRNA Analysis of Tree Soil, Bark, and Trunks

## Overview
A customised bioinformatic pipeline for analyzing 16S rRNA amplicon sequencing data from tree-associated environments, built upon SimpleMetaPipeline and optimised for environmental microbiome data.

## Structure of Research

1. Project Setup & Data Preparation
2. Quality Control & Primer Removal
3. Read Filtering & Trimming
4. Sequence Table Construction
5. LULU Post-Clustering Curation
6. OTU Clustering with Swarm
7. Taxonomic Assignment (IDTAXA (& BLAST))
8. Special Processing for Soil Samples
9. Results Compilation & Export
10. Quality Control & Validation
11. Downstream Analysis Preparation


## Prerequisites

### Software Dependencies
- **Conda** (Miniconda or Anaconda)
- **R** (v4.2.0 or higher) with essential packages: dada2, phyloseq, DECIPHER, lulu
- **Bioinformatics tools**: cutadapt, (BLAST+), swarm
- **Reference databases**: GTDB (primary), SILVA (optional)

### Data Requirements
- Paired-end 16S rRNA sequencing data (FASTQ format)
- Sample metadata CSV file with: sample-id, environment, tree-species, collection-date, location

## Pipeline Execution Steps

### Step 1: Environment Setup and Data Preparation

**Script:** `ControlScriptMethanoProject.r` (initial configuration)

**Actions:**
- Validate FASTQ file integrity and naming conventions
- Parse metadata file and match with sample files
- Create output directory structure
- Initialize logging system


### Step 2: Quality Control and Primer Removal

**Tools:** Cutadapt, FastQC, MultiQC

**Actions:**
- Remove primers and adapters using Cutadapt
- Generate quality reports with FastQC
- Aggregate QC reports with MultiQC
- Remove reads with ambiguous bases (N)
- Trim reads based on quality profiles


### Step 3: Read Filtering and Trimming

**Tool:** DADA2

**Actions:**
- Quality-based filtering and trimming
- Learn error rates from the data
- Dereplication to identify unique sequences
- Sample inference to resolve amplicon sequence variants (ASVs)
- Merge paired-end reads
- Remove chimeric sequences

### Step 4: Sequence Table Construction

**Output:** `data/processed/ASV_table.csv`

**Actions:**
- Create abundance table of ASVs across samples
- Remove singletons and low-abundance sequences
- Apply minimum abundance thresholds:
  - ≥5 reads per sample
  - ≥10 total reads across all samples
  - ≥0.1% relative abundance

### Step 5: LULU Post-Clustering Curation

**Tool:** LULU algorithm

**Actions:**
- First pass: Conservative curation 
- Second pass: Less stringent curation
- Compare curated vs. uncurated tables
- Generate curation statistics and reports

### Step 6: OTU Clustering with Swarm

**Tool:** Swarm

**Actions:**
- Single-linkage clustering with 1bp difference
- Generate OTU table from swarm clusters
- Validate clustering results
- Compare ASV vs. OTU-based approaches

### Step 7: Taxonomic Assignment

**Tools:** IDTAXA (primary), BLAST (fallback)

**Actions:**
- Assign taxonomy using GTDB reference database
- Use BLAST as fallback for unclassified sequences
- Merge results from both methods
- Assign taxonomy at all taxonomic ranks (Phylum to Species)
- Generate assignment statistics and success rates

### Step 8: Special Processing for Soil Samples

**Actions:**
- Apply environment-specific filtering parameters
- Compare soil vs. non-soil processing results
- Generate comparative reports

### Step 9: Results Compilation and Export

**Actions:**
- Run the pipeline
- Compile all processed data into final tables
- Generate summary statistics and visualizations
- Create interactive HTML reports
- Export data in multiple formats for downstream analysis

**Output files:**
- `data/results/ASV_table.csv` - Final abundance table
- `data/results/taxonomy.csv` - Taxonomic assignments
- `data/results/sequence_processing_summary.pdf` - QC report
- `data/results/taxonomic_composition/` - Sample-specific reports


### Step 10: Quality Control and Validation

**Quality checks:**
- Primer removal efficiency reports
- Read retention statistics at each step
- Taxonomic assignment success rates
- Sample-specific processing summaries
- NA value analysis and reporting

### Step 11: Downstream Analysis Preparation

**Output preparation for:**
- Phyloseq object creation
- Diversity analysis (alpha and beta diversity)
- Differential abundance testing
- Visualization and plotting

## Custom Script Execution

After main pipeline completion, run additional analyses:

```bash
# NA value analysis
Rscript scripts/summarize_NA_in_SeqDataTable.R

# Interactive visualization
Rscript scripts/visualize_Methano_SeqDataTable.R

# Taxonomic composition by environment
Rscript scripts/export_taxa_by_sampletype.R

# Curated OTU analysis
Rscript scripts/summarize_NA_curatedOTU.R
```

## Monitoring and Troubleshooting

**Log files:** `logs/pipeline_run_*.log`
**Intermediate files:** `data/processed/` for debugging
**Quality reports:** `data/results/quality_control/`

The pipeline includes comprehensive error handling and can be restarted from any step if interrupted. Each processing step generates validation checks and diagnostic output to ensure data quality throughout the analysis.