#!/usr/bin/env Rscript
# Script 0.3: Read Filtering & Trimming

cat("=== Script 0.3: Read Filtering & Trimming ===\n")

# Load libraries and previous results
library(dada2, quietly = TRUE)

params <- readRDS("data/step_outputs/01_params.rds")
qc_results <- readRDS("data/step_outputs/02_quality_control.rds")

# Quality filtering parameters (adjust these based on quality profiles)
cat("\n=== Quality Filtering Parameters ===\n")
TRUNC_LEN_F <- 240  # Forward read truncation length
TRUNC_LEN_R <- 160  # Reverse read truncation length  
MAX_EE_F <- 2       # Max expected errors for forward reads
MAX_EE_R <- 2       # Max expected errors for reverse reads

cat("Filtering parameters:\n")
cat("  Forward reads: truncLen =", TRUNC_LEN_F, ", maxEE =", MAX_EE_F, "\n")
cat("  Reverse reads: truncLen =", TRUNC_LEN_R, ", maxEE =", MAX_EE_R, "\n")

# Use trimmed files from Step 2
fnFs_trimmed <- qc_results$output_files_F
fnRs_trimmed <- qc_results$output_files_R

# Output paths for filtered files
fnFs_filt <- file.path(params$processed_dir, paste0(params$samples, "_F_filt.fastq.gz"))
fnRs_filt <- file.path(params$processed_dir, paste0(params$samples, "_R_filt.fastq.gz"))
names(fnFs_filt) <- params$samples
names(fnRs_filt) <- params$samples

# Quality filtering and trimming
cat("\n=== Quality Filtering & Trimming ===\n")
cat("Processing", length(params$samples), "samples...\n")

out <- filterAndTrim(fnFs_trimmed, fnFs_filt, fnRs_trimmed, fnRs_filt,
                     truncLen = c(TRUNC_LEN_F, TRUNC_LEN_R),
                     maxN = 0,
                     maxEE = c(MAX_EE_F, MAX_EE_R),
                     truncQ = 2,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = TRUE,
                     verbose = TRUE)

cat("✓ Quality filtering completed\n")
print(out)

# Calculate filtering statistics
reads_in_total <- sum(out[, "reads.in"])
reads_out_total <- sum(out[, "reads.out"])
retention_rate <- round((reads_out_total / reads_in_total) * 100, 1)

cat("\nFiltering Summary:\n")
cat("  Total reads in:", reads_in_total, "\n")
cat("  Total reads out:", reads_out_total, "\n")
cat("  Retention rate:", retention_rate, "%\n")

# Save results
saveRDS(list(
  filtering_results = out,
  filtering_params = list(
    truncLen_F = TRUNC_LEN_F,
    truncLen_R = TRUNC_LEN_R,
    maxEE_F = MAX_EE_F,
    maxEE_R = MAX_EE_R
  ),
  filtered_files_F = fnFs_filt,
  filtered_files_R = fnRs_filt,
  samples_processed = params$samples,
  retention_rate = retention_rate
), "data/step_outputs/03_filtering.rds")

cat("\n✓ Script 0.3 completed successfully!\n")
cat("Filtered files: data/processed/*_filt.fastq.gz\n")
cat("Retention rate:", retention_rate, "%\n")
cat("Next: Rscript scripts/0.4_sequence_table.R\n")
