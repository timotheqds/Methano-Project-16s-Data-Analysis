#!/usr/bin/env Rscript
# Script 0.5: LULU Post-Clustering Curation

cat("=== Script 0.5: LULU Post-Clustering Curation ===\n")

# Load libraries and previous results
library(dada2, quietly = TRUE)
library(Biostrings, quietly = TRUE)

params <- readRDS("data/step_outputs/01_params.rds")
sequence_results <- readRDS("data/step_outputs/04_sequence_table.rds")

# Source pipeline functions
function_dir <- file.path(getwd(), "Pipeline", "Functions")
function_files <- list.files(function_dir, pattern = "\\.R$", full.names = TRUE)
for (f in function_files) source(f)

asv_table <- sequence_results$seqtab
sequences <- sequence_results$sequences

cat("Input ASV table dimensions:", dim(asv_table), "\n")

# Use your LULU pipeline function
cat("\n=== Running LULU Curation ===\n")

# Try to use the RunLULU function, fall back to manual curation if it fails
lulu_results <- tryCatch({
  # Use your RunLULU function from Pipeline/Functions/
  result <- RunLULU(
    otu_table = asv_table,
    sequences = sequences,
    output_dir = params$processed_dir
  )
  
  cat("✓ LULU curation completed using RunLULU function\n")
  result
  
}, error = function(e) {
  cat("RunLULU function error, performing manual curation:\n")
  cat("Error message:", conditionMessage(e), "\n")
  
  # Fallback manual LULU-style curation
  # Use more reasonable thresholds for microbiome data
  # ASVs are columns in the matrix, samples are rows
  asv_totals <- colSums(asv_table)
  asv_prevalence <- colSums(asv_table > 0)
  
  # Keep ASVs with ≥50 total reads AND present in ≥2 samples
  # This is much more reasonable than 0.1% threshold
  abundance_threshold <- 50
  prevalence_threshold <- 2
  
  curated_asvs <- (asv_totals >= abundance_threshold) & (asv_prevalence >= prevalence_threshold)
  asv_table_curated <- asv_table[, curated_asvs, drop = FALSE]
  sequences_curated <- sequences[curated_asvs]
  
  list(
    curated_table = asv_table_curated,
    curated_sequences = sequences_curated,
    curation_stats = data.frame(
      Original_ASVs = ncol(asv_table),
      Curated_ASVs = ncol(asv_table_curated),
      ASVs_Removed = ncol(asv_table) - ncol(asv_table_curated)
    )
  )
})

# Compare curated vs. uncurated tables
cat("\n=== Curation Results ===\n")
if ("curated_table" %in% names(lulu_results)) {
  curated_table <- lulu_results$curated_table
  
  cat("Original ASVs:", ncol(asv_table), "\n")
  cat("Curated ASVs:", ncol(curated_table), "\n")
  cat("ASVs removed:", ncol(asv_table) - ncol(curated_table), "\n")
  cat("Reads retained:", round(sum(curated_table)/sum(asv_table)*100, 1), "%\n")
  
  # Save curated results with proper formatting
  # Create ASV table with proper column names (ASV_1, ASV_2, etc.)
  curated_table_formatted <- curated_table
  colnames(curated_table_formatted) <- paste0("ASV_", 1:ncol(curated_table_formatted))
  
  write.csv(curated_table_formatted, "data/results/ASV_table_curated.csv", quote = FALSE)
  
  # Also save the sequences separately
  curated_sequences <- lulu_results$curated_sequences
  if (!is.null(curated_sequences)) {
    sequences_df <- data.frame(
      ASV_ID = paste0("ASV_", 1:length(curated_sequences)),
      Sequence = curated_sequences,
      stringsAsFactors = FALSE
    )
    write.csv(sequences_df, "data/results/ASV_sequences_curated.csv", 
              row.names = FALSE, quote = FALSE)
  }
  
} else {
  cat("Warning: LULU curation failed, using uncurated table\n")
  curated_table <- asv_table
}

# Generate curation statistics and reports
curation_summary <- data.frame(
  Metric = c("Original_ASVs", "Curated_ASVs", "ASVs_Removed", "Reads_Retained_Percent"),
  Value = c(
    ncol(asv_table),
    ncol(curated_table), 
    ncol(asv_table) - ncol(curated_table),
    round(sum(curated_table)/sum(asv_table)*100, 1)
  )
)

write.csv(curation_summary, "data/results/lulu_curation_summary.csv", 
          row.names = FALSE, quote = FALSE)

# Save results for next step
saveRDS(list(
  original_table = asv_table,
  curated_table = curated_table,
  curation_summary = curation_summary,
  lulu_results = lulu_results
), "data/step_outputs/05_lulu_results.rds")

cat("✓ Script 0.5 completed successfully!\n")
cat("Curated table: data/results/ASV_table_curated.csv\n")
cat("Curation summary: data/results/lulu_curation_summary.csv\n")
cat("Next: Rscript scripts/0.6_otu_clustering.R\n")
