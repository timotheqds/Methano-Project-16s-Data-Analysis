#!/usr/bin/env Rscript

# install_r_packages.R
# Script to install required R packages for MethanoPipeline

packages <- c(
  "dada2",        # Denoising, ASVs
  "phyloseq",     # Microbiome data handling
  "DECIPHER",     # Taxonomic assignment (IDTAXA)
  "lulu",         # Post-clustering curation
  "tidyverse",    # Data manipulation and visualization
  "data.table",   # Fast table handling
  "BiocManager"   # To manage Bioconductor packages
)

# Install missing CRAN packages
cran_pkgs <- packages[!(packages %in% rownames(installed.packages()))]
if (length(cran_pkgs)) {
  install.packages(cran_pkgs, repos = "https://cloud.r-project.org/")
}

# Install Bioconductor dependencies if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

bioc_pkgs <- c("Biostrings", "IRanges", "S4Vectors")  # required by DECIPHER, dada2, phyloseq
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

cat("\nAll required R packages are installed and ready to use.\n")
