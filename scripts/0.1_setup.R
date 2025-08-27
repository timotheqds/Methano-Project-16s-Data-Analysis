#!/usr/bin/env Rscript
# Script 0.1: Project Setup & Data Preparation

# ================================================================================
# MODIFY THESE SETTINGS FOR YOUR DATA:
# ================================================================================

# 1. YOUR SAMPLE NAMES (without _R1/_R2.fastq.gz)
EXPECTED_SAMPLES <- c("A16", "A16B", "A16S", "A19", "A19B", "A19S", 
                     "B6", "B6B", "B6S", "B10", "B10B", "B10S")

# 2. YOUR PRIMER SEQUENCES  
FORWARD_PRIMER <- "AACMGGATTAGATACCCKG"  # 799F
REVERSE_PRIMER <- "ACGTCATCCCCACCTTCC"   # 1193R

# ================================================================================
# PIPELINE CODE 
# ================================================================================

cat("=== Script 0.1: Project Setup & Data Preparation ===\n")

# Load required packages
library(dada2, quietly = TRUE)
library(Biostrings, quietly = TRUE)

# Source pipeline functions
function_dir <- file.path(getwd(), "Pipeline", "Functions")
if (!dir.exists(function_dir)) {
  stop("Pipeline/Functions/ directory not found!")
}

# Source all pipeline functions
function_files <- list.files(function_dir, pattern = "\\.R$", full.names = TRUE)
for (f in function_files) {
  source(f)
}
cat("✓ Loaded", length(function_files), "pipeline functions\n")

# Validate FASTQ file integrity and naming conventions
raw_dir <- file.path(getwd(), "data", "raw")
if (!dir.exists(raw_dir)) {
  stop("data/raw/ directory not found!")
}

# Check if FASTQ files exist and validate naming
missing_files <- c()
for (sample in EXPECTED_SAMPLES) {
  r1_file <- file.path(raw_dir, paste0(sample, "_R1.fastq.gz"))
  r2_file <- file.path(raw_dir, paste0(sample, "_R2.fastq.gz"))
  
  if (!file.exists(r1_file)) missing_files <- c(missing_files, paste0(sample, "_R1.fastq.gz"))
  if (!file.exists(r2_file)) missing_files <- c(missing_files, paste0(sample, "_R2.fastq.gz"))
}

if (length(missing_files) > 0) {
  cat("ERROR: Missing FASTQ files:\n")
  cat(paste(missing_files, collapse = "\n"), "\n")
  stop("Please check your FASTQ files")
}

cat("✓ Found all", length(EXPECTED_SAMPLES), "samples with proper naming\n")

# Create output directory structure
if (!dir.exists("data/step_outputs")) dir.create("data/step_outputs", recursive = TRUE)
if (!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)
if (!dir.exists("data/results")) dir.create("data/results", recursive = TRUE)
if (!dir.exists("logs/steps")) dir.create("logs/steps", recursive = TRUE)

# Initialize logging system
log_file <- file.path("logs", paste0("pipeline_", Sys.Date(), ".log"))
cat("Pipeline started:", as.character(Sys.time()), "\n", file = log_file)

# Parse metadata and validate
metadata_file <- file.path("data", "metadata_template.csv")
if (file.exists(metadata_file)) {
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  cat("✓ Metadata file found\n")
} else {
  cat("Warning: No metadata file found at data/metadata_template.csv\n")
  metadata <- data.frame(Sample = EXPECTED_SAMPLES)
}

# Save configuration parameters
params <- list(
  FWD = FORWARD_PRIMER,
  REV = REVERSE_PRIMER,
  samples = EXPECTED_SAMPLES,
  raw_dir = raw_dir,
  processed_dir = file.path(getwd(), "data", "processed"),
  results_dir = file.path(getwd(), "data", "results"),
  log_file = log_file,
  run_date = Sys.Date()
)

validation <- list(
  samples = EXPECTED_SAMPLES,
  raw_dir = raw_dir,
  all_files_present = TRUE,
  metadata = metadata,
  function_count = length(function_files)
)

saveRDS(params, "data/step_outputs/01_params.rds")
saveRDS(validation, "data/step_outputs/01_validation.rds")

cat("✓ Project setup completed successfully!\n")
cat("✓ Output directories created\n")
cat("✓ Logging system initialized\n")
cat("Next: Rscript scripts/0.2_quality_control.R\n")
