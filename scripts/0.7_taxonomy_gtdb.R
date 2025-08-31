#!/usr/bin/env Rscript
# Script 0.7: Taxonomic Assignment (GTDB + SILVA fallback)

cat("=== Script 0.7: Taxonomic Assignment with GTDB ===\n")

# Load libraries and functions
library(dada2, quietly = TRUE)
library(Biostrings, quietly = TRUE)
library(DECIPHER, quietly = TRUE)

# Source pipeline functions
function_dir <- file.path(getwd(), "Pipeline", "Functions")
function_files <- list.files(function_dir, pattern = "\\.R$", full.names = TRUE)
for (f in function_files) source(f)

params <- readRDS("data/step_outputs/01_params.rds")
otu_results <- readRDS("data/step_outputs/06_otu_clustering.rds")
sequence_results <- readRDS("data/step_outputs/04_sequence_table.rds")

# Get sequences for taxonomy assignment
sequences <- sequence_results$sequences
otu_table <- otu_results$otu_table

cat("Assigning taxonomy for", length(sequences), "sequences\n")

# Primary assignment with GTDB
cat("\n=== Primary Taxonomic Assignment with GTDB ===\n")

# Check if GTDB database is available
gtdb_db_path <- "reference_databases/GTDB_r226-mod_April2025.RData"

taxonomy_gtdb <- NULL

if (file.exists(gtdb_db_path)) {
  cat("Using GTDB r226 database for taxonomy assignment\n")
  
  taxonomy_gtdb <- tryCatch({
    # Load GTDB training set
    cat("Loading GTDB database...\n")
    load(gtdb_db_path)  # This should load the trainingSet object
    
    # Convert sequences to DNAStringSet if needed
    if (!inherits(sequences, "DNAStringSet")) {
      sequences_dna <- DNAStringSet(sequences)
    } else {
      sequences_dna <- sequences
    }
    
    cat("Running IDTAXA classification with GTDB...\n")
    # Use IDTAXA for classification with GTDB (lower threshold for better success)
    ids <- IdTaxa(sequences_dna, trainingSet, strand="both", 
                  threshold=30, processors=NULL)
    
    # Convert IDTAXA results to data frame
    taxonomy_df <- data.frame(
      Kingdom = character(length(ids)),
      Phylum = character(length(ids)),
      Class = character(length(ids)),
      Order = character(length(ids)),
      Family = character(length(ids)),
      Genus = character(length(ids)),
      Species = character(length(ids)),
      Method = "GTDB",
      Confidence = numeric(length(ids)),
      stringsAsFactors = FALSE
    )
    
    # Parse IDTAXA results
    for (i in seq_along(ids)) {
      if (length(ids[[i]]$taxon) > 0) {
        # Extract taxonomic levels and confidences
        taxa_levels <- ids[[i]]$taxon
        confidences <- ids[[i]]$confidence
        
        # GTDB hierarchy: Root | Domain | Phylum | Class | Order | Family | Genus | Species
        # Map to standard levels
        for (j in seq_along(taxa_levels)) {
          taxon <- taxa_levels[j]
          conf <- confidences[j]
          
          # Skip Root level
          if (taxon == "Root") next
          
          # Map based on position in hierarchy (after Root)
          level_pos <- j - 1  # Subtract 1 to skip Root
          
          if (level_pos == 1) {  # Domain (Kingdom)
            taxonomy_df$Kingdom[i] <- taxon
            taxonomy_df$Confidence[i] <- max(taxonomy_df$Confidence[i], conf)
          } else if (level_pos == 2) {  # Phylum
            taxonomy_df$Phylum[i] <- taxon
            taxonomy_df$Confidence[i] <- max(taxonomy_df$Confidence[i], conf)
          } else if (level_pos == 3) {  # Class
            taxonomy_df$Class[i] <- taxon
            taxonomy_df$Confidence[i] <- max(taxonomy_df$Confidence[i], conf)
          } else if (level_pos == 4) {  # Order
            taxonomy_df$Order[i] <- taxon
            taxonomy_df$Confidence[i] <- max(taxonomy_df$Confidence[i], conf)
          } else if (level_pos == 5) {  # Family
            taxonomy_df$Family[i] <- taxon
            taxonomy_df$Confidence[i] <- max(taxonomy_df$Confidence[i], conf)
          } else if (level_pos == 6) {  # Genus
            taxonomy_df$Genus[i] <- taxon
            taxonomy_df$Confidence[i] <- max(taxonomy_df$Confidence[i], conf)
          } else if (level_pos == 7) {  # Species
            taxonomy_df$Species[i] <- taxon
            taxonomy_df$Confidence[i] <- max(taxonomy_df$Confidence[i], conf)
          }
        }
      }
    }
    
    # Clean up empty entries
    taxonomy_df[taxonomy_df == ""] <- NA
    
    # Identify chloroplasts and mitochondria
    chloroplast_idx <- grepl("Chloroplast", taxonomy_df$Class, ignore.case = TRUE) | 
                      grepl("Chloroplast", taxonomy_df$Order, ignore.case = TRUE)
    mitochondria_idx <- grepl("Mitochondria", taxonomy_df$Family, ignore.case = TRUE) | 
                       grepl("Mitochondria", taxonomy_df$Genus, ignore.case = TRUE)
    
    # Flag these sequences
    taxonomy_df$Method[chloroplast_idx] <- "GTDB_Chloroplast"
    taxonomy_df$Method[mitochondria_idx] <- "GTDB_Mitochondria"
    
    n_chloroplasts <- sum(chloroplast_idx, na.rm = TRUE)
    n_mitochondria <- sum(mitochondria_idx, na.rm = TRUE)
    n_bacterial <- nrow(taxonomy_df) - n_chloroplasts - n_mitochondria
    
    cat("✓ GTDB taxonomy assignment completed\n")
    cat("Assigned taxonomy for", nrow(taxonomy_df), "sequences\n")
    cat("  - Bacterial sequences:", n_bacterial, "\n")
    cat("  - Chloroplast sequences:", n_chloroplasts, "\n") 
    cat("  - Mitochondrial sequences:", n_mitochondria, "\n")
    
    # Summary of assignments
    kingdom_counts <- table(taxonomy_df$Kingdom, useNA = "ifany")
    cat("Kingdom distribution:\n")
    print(kingdom_counts)
    
    phylum_counts <- table(taxonomy_df$Phylum, useNA = "ifany")
    cat("\nTop 10 Phyla:\n")
    print(head(sort(phylum_counts, decreasing = TRUE), 10))
    
    taxonomy_df
    
  }, error = function(e) {
    cat("GTDB taxonomy assignment failed:", e$message, "\n")
    cat("Error details:", toString(e), "\n")
    NULL
  })
  
} else {
  cat("GTDB database not found at:", gtdb_db_path, "\n")
  cat("Please ensure the GTDB database is properly installed\n")
}

# Fallback to SILVA for unclassified sequences
cat("\n=== Fallback Taxonomic Assignment with SILVA ===\n")

taxonomy_final <- taxonomy_gtdb

if (!is.null(taxonomy_gtdb)) {
  # Check which bacterial sequences need additional classification (exclude chloroplasts/mitochondria)
  bacterial_seqs <- !grepl("Chloroplast|Mitochondria", taxonomy_gtdb$Method)
  unclassified <- bacterial_seqs & (is.na(taxonomy_gtdb$Genus) | taxonomy_gtdb$Genus == "")
  n_unclassified <- sum(unclassified, na.rm = TRUE)
  
  if (n_unclassified > 0) {
    cat("Found", n_unclassified, "bacterial sequences unclassified by GTDB\n")
    cat("Attempting SILVA classification for these sequences\n")
    
    # Try SILVA for unclassified sequences
    silva_db_path <- "reference_databases/silva_nr99_v138.1_train_set.fa.gz"
    silva_species_path <- "reference_databases/silva_species_assignment_v138.1.fa.gz"
    
    if (file.exists(silva_db_path)) {
      silva_results <- tryCatch({
        # Get unclassified sequences
        unclass_sequences <- sequences[unclassified]
        
        # Assign taxonomy using SILVA
        taxa_silva <- assignTaxonomy(unclass_sequences, silva_db_path, multithread = TRUE)
        
        # Add species assignment if available
        if (file.exists(silva_species_path)) {
          taxa_silva <- addSpecies(taxa_silva, silva_species_path)
        }
        
        # Update the final taxonomy with SILVA results
        taxonomy_final[unclassified, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")] <- 
          taxa_silva[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")]
        
        if ("Species" %in% colnames(taxa_silva)) {
          taxonomy_final$Species[unclassified] <- taxa_silva[, "Species"]
        }
        
        taxonomy_final$Method[unclassified] <- "SILVA_fallback"
        taxonomy_final$Confidence[unclassified] <- 0.6  # Lower confidence for fallback
        
        cat("✓ SILVA fallback completed for", sum(!is.na(taxa_silva[, "Genus"])), "sequences\n")
        
      }, error = function(e) {
        cat("SILVA fallback failed:", e$message, "\n")
      })
    }
  }
}

# Final fallback for still unclassified
if (!is.null(taxonomy_final)) {
  still_unclassified <- is.na(taxonomy_final$Genus) | taxonomy_final$Genus == ""
  n_still_unclassified <- sum(still_unclassified, na.rm = TRUE)
  
  if (n_still_unclassified > 0) {
    cat("Applying basic classification to", n_still_unclassified, "remaining sequences\n")
    
    taxonomy_final$Kingdom[still_unclassified] <- "Bacteria"
    taxonomy_final$Phylum[still_unclassified] <- "Unknown_Bacteria"
    taxonomy_final$Class[still_unclassified] <- "Unknown_Bacteria"
    taxonomy_final$Order[still_unclassified] <- "Unknown_Bacteria"
    taxonomy_final$Family[still_unclassified] <- "Unknown_Bacteria"
    taxonomy_final$Genus[still_unclassified] <- "Unknown_Bacteria"
    taxonomy_final$Method[still_unclassified] <- "Fallback"
    taxonomy_final$Confidence[still_unclassified] <- 0.1
  }
}

# Create final taxonomy results
if (!is.null(taxonomy_final)) {
  # Summary statistics
  cat("\n=== Final Taxonomy Summary ===\n")
  cat("Total sequences:", nrow(taxonomy_final), "\n")
  
  method_counts <- table(taxonomy_final$Method)
  cat("Assignment methods:\n")
  print(method_counts)
  
  # Count classified vs unclassified at different levels (excluding chloroplasts/mitochondria)
  bacterial_only <- !grepl("Chloroplast|Mitochondria", taxonomy_final$Method)
  bacterial_count <- sum(bacterial_only, na.rm = TRUE)
  
  if (bacterial_count > 0) {
    classified_genus <- sum(bacterial_only & !is.na(taxonomy_final$Genus) & taxonomy_final$Genus != "Unknown_Bacteria", na.rm = TRUE)
    classified_family <- sum(bacterial_only & !is.na(taxonomy_final$Family) & taxonomy_final$Family != "Unknown_Bacteria", na.rm = TRUE)
    classified_phylum <- sum(bacterial_only & !is.na(taxonomy_final$Phylum) & taxonomy_final$Phylum != "Unknown_Bacteria", na.rm = TRUE)
    
    cat("\nBacterial classification success rates (excluding chloroplasts/mitochondria):\n")
    cat("Total bacterial sequences:", bacterial_count, "\n")
    cat("Phylum level:", round(classified_phylum/bacterial_count*100, 1), "%\n")
    cat("Family level:", round(classified_family/bacterial_count*100, 1), "%\n")
    cat("Genus level:", round(classified_genus/bacterial_count*100, 1), "%\n")
  }
  
  # Create results object
  taxonomy_results <- list(
    taxonomy_table = taxonomy_final,
    otu_table = otu_table,
    sequences = sequences,
    method_summary = method_counts,
    classification_rates = c(
      phylum = classified_phylum/nrow(taxonomy_final),
      family = classified_family/nrow(taxonomy_final),
      genus = classified_genus/nrow(taxonomy_final)
    )
  )
  
  # Save results
  output_file <- "data/step_outputs/07_taxonomy_gtdb.rds"
  saveRDS(taxonomy_results, output_file)
  cat("✓ Results saved to:", output_file, "\n")
  
  # Also save a CSV for easy viewing
  csv_file <- "data/results/taxonomy_assignments_gtdb.csv"
  if (!dir.exists("data/results")) dir.create("data/results", recursive = TRUE)
  
  # We need to map taxonomy back to OTU representatives
  # The OTU table has 215 rows but we classified 3449 sequences
  # We need to extract taxonomy for sequences that correspond to OTU representatives
  
  if (nrow(otu_table) == length(sequences)) {
    # Direct mapping - OTU table matches sequence count
    final_table <- cbind(otu_table, taxonomy_final[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Method", "Confidence")])
  } else {
    # We need to map sequences to OTU representatives
    # For now, save just the taxonomy table
    cat("Note: OTU table (", nrow(otu_table), " rows) doesn't match taxonomy table (", nrow(taxonomy_final), " rows)\n")
    cat("Saving taxonomy table separately\n")
    final_table <- taxonomy_final
  }
  
  write.csv(final_table, csv_file, row.names = TRUE)
  cat("✓ Taxonomy table saved to:", csv_file, "\n")
  
  cat("\n=== Taxonomy Assignment with GTDB Complete ===\n")
  
} else {
  cat("❌ Taxonomy assignment failed completely\n")
  stop("No taxonomy assignment method succeeded")
}
