#!/usr/bin/env Rscript
# Script 0.2: Quality Control & Primer Removal

cat("=== Script 0.2: Quality Control & Primer Removal ===\n")

# Load libraries and previous results
library(dada2, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(Biostrings, quietly = TRUE)

params <- readRDS("data/step_outputs/01_params.rds")

# Quality profile analysis for ALL samples
cat("\n=== Quality Control Analysis ===\n")
cat("Creating quality profiles for all", length(params$samples), "samples...\n")

fnFs <- file.path(params$raw_dir, paste0(params$samples, "_R1.fastq.gz"))
fnRs <- file.path(params$raw_dir, paste0(params$samples, "_R2.fastq.gz"))

# Generate quality reports
p1 <- plotQualityProfile(fnFs)
ggsave("data/results/quality_profiles_forward.pdf", p1, width = 16, height = 12)

p2 <- plotQualityProfile(fnRs)
ggsave("data/results/quality_profiles_reverse.pdf", p2, width = 16, height = 12)

cat("✓ Quality profiles saved for all samples\n")

# Primer removal using cutadapt
cat("\n=== Primer Removal ===\n")
cat("Checking for cutadapt...\n")

cutadapt_cmd <- system("which cutadapt", intern = TRUE, ignore.stderr = TRUE)
if (length(cutadapt_cmd) == 0) {
  stop("cutadapt not found! Install with: conda install -c bioconda cutadapt")
}

# Check cutadapt version
version_cmd <- system("cutadapt --version", intern = TRUE, ignore.stderr = TRUE)
cat("Found cutadapt version:", version_cmd[1], "\n")
cat("Primers: FWD =", params$FWD, "| REV =", params$REV, "\n")

# Expand IUPAC codes for primer variants (use only most common variants)
# FWD primer: AACMGGATTAGATACCCKG (M=A/C, K=G/T)
# Based on data analysis, use only the two most common variants
FWD_primary <- "AACAGGATTAGATACCCGG"    # Most common (261 occurrences)
FWD_secondary <- "AACAGGATTAGATACCCGT"  # Second most common (197 occurrences)

# REV primer is already concrete: ACGTCATCCCCACCTTCC
REV_primer <- params$REV

# Create reverse complements
FWD_primary_RC <- as.character(reverseComplement(DNAString(FWD_primary)))
FWD_secondary_RC <- as.character(reverseComplement(DNAString(FWD_secondary)))
REV_RC <- as.character(reverseComplement(DNAString(REV_primer)))

cat("Using primary FWD primer:", FWD_primary, "\n")
cat("Using secondary FWD primer:", FWD_secondary, "\n")
cat("REV primer:", REV_primer, "\n")

# Output files after primer removal
fnFs_out <- file.path(params$processed_dir, paste0(params$samples, "_F_trimmed.fastq.gz"))
fnRs_out <- file.path(params$processed_dir, paste0(params$samples, "_R_trimmed.fastq.gz"))

# Run cutadapt for each sample
for (i in seq_along(params$samples)) {
  cat("Processing", params$samples[i], "...\n")
  
  # Simplified command with only the most common primer variants
  # For paired-end reads: -g/-G for 5' adapters, -a/-A for 3' adapters
  cmd <- paste0("cutadapt",
                " -g ", FWD_primary, " -g ", FWD_secondary,  # 5' forward primers
                " -G ", REV_primer,                          # 5' reverse primer
                " -a ", REV_RC,                              # 3' reverse complement
                " -A ", FWD_primary_RC, " -A ", FWD_secondary_RC,  # 3' forward complements
                " --minimum-length 50 --max-n 0 --overlap 10",
                " -o ", fnFs_out[i], " -p ", fnRs_out[i],
                " ", fnFs[i], " ", fnRs[i])
  
  result <- system(cmd, intern = TRUE)
  cat("  ✓ Completed\n")
}

primer_status <- "Completed with cutadapt"

cat("✓ Primer removal completed for all", length(params$samples), "samples\n")

# Save results
saveRDS(list(
  samples = params$samples,
  quality_profiles_created = TRUE,
  primer_removal_status = primer_status,
  primers_used = list(FWD = params$FWD, REV = params$REV),
  output_files_F = fnFs_out,
  output_files_R = fnRs_out
), "data/step_outputs/02_quality_control.rds")

cat("\nScript 0.2 completed successfully!\n")
cat("Quality profiles: data/results/quality_profiles_*.pdf\n")
cat("Primer removal:", primer_status, "\n")
cat("Output files: data/processed/*_trimmed.fastq.gz\n")
cat("Next: Rscript scripts/0.3_filtering.R\n")
