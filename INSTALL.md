# MethanoPipeline: Prerequisites and Installation Guide

## Quick Start for New Users

### Prerequisites
Before cloning this repository, install:

1. **R** (v4.2.0 or higher): https://cran.r-project.org/
2. **Conda** (Miniconda recommended): https://docs.conda.io/en/latest/miniconda.html
3. **Git**: https://git-scm.com/downloads

### Installation Steps

```bash
# 1. Clone the repository
git clone https://github.com/timotheqds/Methano-Project-16s-Data-Analysis.git
cd Methano-Project-16s-Data-Analysis

# 2. Set up the conda environment (includes all bioinformatics tools)
conda env create -f environment.yml
conda activate methanopipeline

# 3. Install R packages
Rscript scripts/install_r_packages.R

# 4. Download reference databases (optional - can be skipped for testing)
Rscript scripts/download_databases.R

# 5. Verify installation
Rscript scripts/verify_installation.R
```

### Data Setup

```bash
# Place your FASTQ files in data/raw/
# Files should be named: SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz

# Create your metadata file in data/metadata_template.csv
# See data/metadata_template.csv for format
```

### Running the Pipeline

```bash
# Run all steps sequentially
Rscript scripts/step_controller.R run step1
Rscript scripts/step_controller.R run step2
Rscript scripts/step_controller.R run step3
Rscript scripts/step_controller.R run step4

# Or run individual steps for testing
Rscript scripts/steps/step1_setup_validation.R
```

### Expected Output Structure
```
data/
├── raw/                    # Your FASTQ files
├── processed/              # Intermediate files
├── results/                # Final results
└── step_outputs/           # Step-by-step outputs

logs/
└── steps/                  # Execution logs
```

### Customization for Your Data

Edit these files for your specific dataset:
- `data/metadata_template.csv` - Your sample information
- `scripts/steps/step1_setup_validation.R` - Sample names and primer sequences
- See README.md for detailed customization options

### Troubleshooting

Common issues and solutions are documented in `docs/troubleshooting.md`
