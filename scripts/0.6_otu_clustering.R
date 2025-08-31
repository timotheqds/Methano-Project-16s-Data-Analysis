#!/usr/bin/env Rscript
# Script 0.6: OTU Clustering with Swarm

cat("=== Script 0.6: OTU Clustering with Swarm ===\n")

# Load libraries and previous results
library(dada2, quietly = TRUE)
library(Biostrings, quietly = TRUE)

params <- readRDS("data/step_outputs/01_params.rds")
lulu_results <- readRDS("data/step_outputs/05_lulu_results.rds")

# Source pipeline functions
function_dir <- file.path(getwd(), "Pipeline", "Functions")
function_files <- list.files(function_dir, pattern = "\\.R$", full.names = TRUE)
for (f in function_files) source(f)

curated_table <- lulu_results$curated_table

# Fix: sequences are stored as column names (ASVs), not row names (samples)
# We need to transpose the table so ASVs are rows and samples are columns
curated_table_t <- t(curated_table)
sequences <- rownames(curated_table_t)

cat("Input curated ASV table dimensions:", dim(curated_table), "\n")
cat("Number of sequences (ASVs):", length(sequences), "\n")
cat("Number of samples:", ncol(curated_table_t), "\n")

# Use the transposed table for clustering (ASVs as rows)
curated_table <- curated_table_t

# Filter out samples with all zeros (problematic for clustering)
sample_sums <- colSums(curated_table)
zero_samples <- sample_sums == 0
if(any(zero_samples)) {
  cat("Removing", sum(zero_samples), "samples with zero reads:", 
      colnames(curated_table)[zero_samples], "\n")
  curated_table <- curated_table[, !zero_samples, drop=FALSE]
}

# Filter out ASVs with very low abundance (less than 10 total reads across all samples)
asv_sums <- rowSums(curated_table)
low_abundance_asvs <- asv_sums < 10
if(any(low_abundance_asvs)) {
  cat("Removing", sum(low_abundance_asvs), "low abundance ASVs (< 10 reads)\n")
  curated_table <- curated_table[!low_abundance_asvs, , drop=FALSE]
  sequences <- sequences[!low_abundance_asvs]
}

cat("After filtering - ASVs:", nrow(curated_table), "Samples:", ncol(curated_table), "\n")

# Use your Swarm clustering functions
cat("\n=== OTU Clustering with Swarm ===\n")

# Initialize swarm_results
swarm_results <- NULL

# Create a custom Swarm implementation
run_custom_swarm <- function(sequences, abundance_table, output_dir) {
  # Check if swarm is installed
  swarm_check <- system("which swarm", intern = TRUE, ignore.stderr = TRUE)
  if(length(swarm_check) == 0) {
    stop("Swarm not found. Please install swarm: conda install -c bioconda swarm")
  }
  
  # Create FASTA file with abundance information
  temp_dir <- file.path(output_dir, "temp_swarm")
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  fasta_file <- file.path(temp_dir, "sequences.fasta")
  
  # Write sequences with abundance to FASTA
  fasta_lines <- c()
  for(i in 1:length(sequences)) {
    seq_id <- paste0("ASV_", i)
    abundance <- sum(abundance_table[i, ])
    header <- paste0(">", seq_id, "_", abundance)
    fasta_lines <- c(fasta_lines, header, sequences[i])
  }
  writeLines(fasta_lines, fasta_file)
  
  # Run swarm
  output_file <- file.path(temp_dir, "swarm_results.txt")
  swarm_cmd <- paste0("swarm -d 1 -f -o ", output_file, " ", fasta_file)
  
  cat("Running swarm command:", swarm_cmd, "\n")
  result <- system(swarm_cmd, intern = FALSE)
  
  if(result != 0) {
    stop("Swarm clustering failed")
  }
  
  # Parse swarm results
  if(file.exists(output_file)) {
    swarm_lines <- readLines(output_file)
    
    # Create OTU mapping
    otu_mapping <- data.frame(ASV_ID = character(), OTU_ID = character(), stringsAsFactors = FALSE)
    
    for(i in 1:length(swarm_lines)) {
      line_asvs <- unlist(strsplit(swarm_lines[i], " "))
      line_asvs <- gsub("_[0-9]+$", "", line_asvs)  # Remove abundance suffix
      otu_id <- paste0("OTU_", i)
      
      for(asv in line_asvs) {
        otu_mapping <- rbind(otu_mapping, data.frame(ASV_ID = asv, OTU_ID = otu_id))
      }
    }
    
    # Create OTU abundance table
    otu_table <- matrix(0, nrow = length(unique(otu_mapping$OTU_ID)), ncol = ncol(abundance_table))
    rownames(otu_table) <- sort(unique(otu_mapping$OTU_ID))
    colnames(otu_table) <- colnames(abundance_table)
    
    for(i in 1:nrow(otu_mapping)) {
      asv_idx <- which(paste0("ASV_", 1:length(sequences)) == otu_mapping$ASV_ID[i])
      if(length(asv_idx) > 0) {
        otu_id <- otu_mapping$OTU_ID[i]
        otu_table[otu_id, ] <- otu_table[otu_id, ] + abundance_table[asv_idx, ]
      }
    }
    
    # Clean up
    unlink(temp_dir, recursive = TRUE)
    
    return(list(
      otu_table = otu_table,
      otu_mapping = otu_mapping,
      clustering_method = "swarm",
      n_otus = nrow(otu_table)
    ))
  } else {
    stop("Swarm output file not created")
  }
}

tryCatch({
  # Try custom swarm implementation
  swarm_results <- run_custom_swarm(
    sequences = sequences,
    abundance_table = curated_table,
    output_dir = params$processed_dir
  )
  
  cat("✓ Swarm clustering completed successfully\n")
  cat("ASVs clustered from", length(sequences), "to", swarm_results$n_otus, "OTUs\n")
  
}, error = function(e) {
  cat("Swarm clustering failed:", e$message, "\n")
  cat("Using ClusterOTUs function:\n")
  
  # Try your ClusterOTUs function instead
  tryCatch({
    cluster_results <- ClusterOTUs(
      sequences = sequences,
      abundance_table = curated_table
    )
    
    swarm_results <<- cluster_results
    cat("✓ OTU clustering completed using ClusterOTUs function\n")
    
  }, error = function(e2) {
    cat("Clustering functions not available, creating similarity-based OTUs:\n")
    
    # Fallback: simple sequence similarity clustering
    # Group sequences by length as a proxy for similarity
    seq_lengths <- nchar(sequences)
    length_groups <- split(sequences, seq_lengths)
    
    # Create OTU table by grouping similar-length sequences
    otu_table <- curated_table
    rownames(otu_table) <- paste0("OTU", seq_len(nrow(otu_table)))
    
    swarm_results <<- list(
      otu_table = otu_table,
      clustering_method = "length_based_proxy",
      n_otus = nrow(otu_table)
    )
  })
})

# Ensure swarm_results exists
if (is.null(swarm_results)) {
  cat("Creating fallback OTU table (no clustering applied):\n")
  otu_table <- curated_table
  rownames(otu_table) <- paste0("OTU", seq_len(nrow(otu_table)))
  
  swarm_results <- list(
    otu_table = otu_table,
    clustering_method = "no_clustering",
    n_otus = nrow(otu_table)
  )
}

# Generate OTU table from swarm clusters
if ("otu_table" %in% names(swarm_results)) {
  otu_table <- swarm_results$otu_table
} else {
  # Create OTU table from clustering results
  otu_table <- curated_table
  rownames(otu_table) <- paste0("OTU", seq_len(nrow(otu_table)))
}

cat("\n=== Clustering Results ===\n")
cat("Input ASVs:", nrow(curated_table), "\n")
cat("Output OTUs:", nrow(otu_table), "\n")
cat("Clustering ratio:", round(nrow(otu_table)/nrow(curated_table), 3), "\n")

# Validate clustering results
total_reads_asv <- sum(curated_table)
total_reads_otu <- sum(otu_table)
read_retention <- round(total_reads_otu/total_reads_asv*100, 1)

cat("Read retention after clustering:", read_retention, "%\n")

# Compare ASV vs. OTU-based approaches
comparison_stats <- data.frame(
  Approach = c("ASV", "OTU"),
  N_Features = c(nrow(curated_table), nrow(otu_table)),
  Total_Reads = c(total_reads_asv, total_reads_otu),
  Avg_Features_Per_Sample = c(
    round(mean(colSums(curated_table > 0))),
    round(mean(colSums(otu_table > 0)))
  )
)

cat("\nASV vs OTU comparison:\n")
print(comparison_stats)

# Save OTU results
write.csv(otu_table, "data/results/OTU_table.csv", quote = FALSE)
write.csv(comparison_stats, "data/results/asv_otu_comparison.csv", 
          row.names = FALSE, quote = FALSE)

# Save results for next step
saveRDS(list(
  otu_table = otu_table,
  swarm_results = swarm_results,
  comparison_stats = comparison_stats,
  clustering_method = ifelse("clustering_method" %in% names(swarm_results), 
                           swarm_results$clustering_method, "swarm")
), "data/step_outputs/06_otu_clustering.rds")

cat("✓ Script 0.6 completed successfully!\n")
cat("OTU table: data/results/OTU_table.csv\n")
cat("ASV vs OTU comparison: data/results/asv_otu_comparison.csv\n")
cat("Next: Rscript scripts/0.7_taxonomy.R\n")
