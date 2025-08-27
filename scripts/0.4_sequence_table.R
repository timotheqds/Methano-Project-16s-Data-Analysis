#!/usr/bin/env Rscript
# Script 0.4: Sequence Table Construction

cat("=== Script 0.4: Sequence Table Construction ===\n")

# Load libraries and previous results
library(dada2, quietly = TRUE)

params <- readRDS("data/step_outputs/01_params.rds")
filtering_results <- readRDS("data/step_outputs/03_filtering.rds")

# Get filtered file paths
fnFs_filt <- filtering_results$filtered_files_F
fnRs_filt <- filtering_results$filtered_files_R

# Check if all filtered files exist
if (!all(file.exists(fnFs_filt)) || !all(file.exists(fnRs_filt))) {
  stop("Some filtered files are missing!")
}

cat("✓ All", length(fnFs_filt), "filtered files found\n")

# Learn error rates
cat("\n=== Learning Error Rates ===\n")
cat("Forward reads error learning...\n")
errF <- learnErrors(fnFs_filt, multithread = TRUE, verbose = 1)

cat("Reverse reads error learning...\n") 
errR <- learnErrors(fnRs_filt, multithread = TRUE, verbose = 1)

cat("✓ Error learning completed\n")

# Sample inference (denoising)
cat("\n=== Sample Inference (Denoising) ===\n")
dadaFs <- dada(fnFs_filt, err = errF, multithread = TRUE, verbose = 1)
dadaRs <- dada(fnRs_filt, err = errR, multithread = TRUE, verbose = 1)

cat("✓ Denoising completed\n")

# Merge paired reads  
cat("\n=== Merging Paired Reads ===\n")
mergers <- mergePairs(dadaFs, fnFs_filt, dadaRs, fnRs_filt, verbose = TRUE)

cat("✓ Read merging completed\n")

# Construct sequence table
cat("\n=== Constructing Sequence Table ===\n")
seqtab <- makeSequenceTable(mergers)
cat("Sequence table dimensions:", dim(seqtab), "\n")
cat("  Samples:", nrow(seqtab), "\n")
cat("  ASVs:", ncol(seqtab), "\n")

# Inspect distribution of sequence lengths
seq_lengths <- nchar(getSequences(seqtab))
cat("Sequence length distribution:\n")
print(table(seq_lengths))

# Remove chimeras
cat("\n=== Removing Chimeras ===\n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
cat("Chimera removal summary:\n")
cat("  Before chimera removal:", ncol(seqtab), "ASVs\n")
cat("  After chimera removal:", ncol(seqtab_nochim), "ASVs\n")
cat("  Reads retained:", round(sum(seqtab_nochim)/sum(seqtab)*100, 1), "%\n")

# Track reads through the pipeline
cat("\n=== Read Tracking Summary ===\n")
getN <- function(x) sum(getUniques(x))

# Create tracking table
track <- cbind(
  filtering_results$filtering_results[, "reads.out"],
  sapply(dadaFs, getN),
  sapply(dadaRs, getN), 
  sapply(mergers, getN),
  rowSums(seqtab_nochim)
)

colnames(track) <- c("filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- params$samples

print(track)

# Save results for next step (LULU curation)
saveRDS(list(
  seqtab = seqtab_nochim,
  sequences = getSequences(seqtab_nochim),
  track = track,
  error_rates = list(errF = errF, errR = errR),
  dada_results = list(dadaFs = dadaFs, dadaRs = dadaRs),
  mergers = mergers,
  samples = params$samples
), "data/step_outputs/04_sequence_table.rds")

cat("\n✓ Script 0.4 completed successfully!\n")
cat("Sequence table construction completed:\n")
cat("  Final ASVs:", ncol(seqtab_nochim), "\n")
cat("  Total reads retained:", sum(seqtab_nochim), "\n")
cat("  Average reads per sample:", round(mean(rowSums(seqtab_nochim))), "\n")
cat("Next: Rscript scripts/0.5_lulu_curation.R\n")
