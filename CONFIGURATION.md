# MethanoPipeline Configuration Template

## Step-by-Step Customization Guide

### STEP 1: Update Sample Names

**File:** `scripts/steps/step1_setup_validation.R`  
**Line:** ~31

```r
# MODIFY THIS SECTION FOR YOUR DATA:
expected_samples <- c("SAMPLE1", "SAMPLE2", "SAMPLE3")  # Replace with your sample names
```

**Requirements:**
- Your FASTQ files must be named: `SAMPLENAME_R1.fastq.gz` and `SAMPLENAME_R2.fastq.gz`
- Place all FASTQ files in `data/raw/` directory

### STEP 2: Update Primer Sequences

**File:** `scripts/steps/step1_setup_validation.R`  
**Line:** ~86

```r
# MODIFY THESE PRIMERS FOR YOUR DATA:
FWD <- "AACMGGATTAGATACCCKG"  # Replace with your forward primer
REV <- "ACGTCATCCCCACCTTCC"   # Replace with your reverse primer
```

**Common primer pairs:**
- V3-V4: FWD="CCTACGGGNGGCWGCAG", REV="GACTACHVGGGTATCTAATCC"
- V4: FWD="GTGCCAGCMGCCGCGGTAA", REV="GGACTACHVGGGTWTCTAAT"

### STEP 3: Update DADA2 Parameters

**File:** `scripts/steps/step1_setup_validation.R`  
**Line:** ~95

```r
# MODIFY THESE PARAMETERS BASED ON YOUR QUALITY PROFILES:
params <- list(
  truncLen = c(240, 160),               # [R1_length, R2_length] - ADJUST AFTER VIEWING QUALITY PLOTS
  trimLeft = c(0, 0),                   # Bases to trim from start [R1, R2]
  maxEE = c(2, 2),                      # Max expected errors [R1, R2]
  DesiredSequenceLengthRange = c(250, 450)  # Expected amplicon size range
)
```

**Parameter Guidelines:**
1. **Run Step 2 first** to see quality profiles
2. Set `truncLen` where quality drops below Q20
3. Adjust `DesiredSequenceLengthRange` based on your primers
4. Start with `maxEE = c(2,2)`, increase if too few reads pass

### STEP 4: Create Metadata File

**File:** `data/metadata_template.csv`

**Required format:**
```csv
sample_id,sample_type,group,description
SAMPLE1,Soil,Site1,Description of sample 1
SAMPLE2,Water,Site1,Description of sample 2
SAMPLE3,Sediment,Site2,Description of sample 3
```

**Requirements:**
- `sample_id` must match your FASTQ file names exactly
- Add any additional columns you need for your analysis

### STEP 5: Test Configuration

```bash
# Test your configuration
Rscript scripts/steps/step1_setup_validation.R

# If successful, you should see:
# ✓ All expected paired-end files present
# ✓ Metadata template found
# ✓ Step 1 completed successfully!
```

## Quick Configuration Examples

### Example 1: Soil Microbiome Study
```r
# Sample names
expected_samples <- c("Soil_A1", "Soil_A2", "Soil_B1", "Soil_B2")

# V4 primers
FWD <- "GTGCCAGCMGCCGCGGTAA"  # 515F
REV <- "GGACTACHVGGGTWTCTAAT"  # 806R

# Parameters for V4 region
truncLen = c(250, 200)
DesiredSequenceLengthRange = c(240, 260)
```

### Example 2: Marine Samples
```r
# Sample names  
expected_samples <- c("Marine_1", "Marine_2", "Marine_3")

# V3-V4 primers
FWD <- "CCTACGGGNGGCWGCAG"    # 341F
REV <- "GACTACHVGGGTATCTAATCC"  # 785R

# Parameters for V3-V4 region
truncLen = c(280, 200)
DesiredSequenceLengthRange = c(400, 460)
```

## Validation Checklist

Before running the full pipeline:

- [ ] Updated `expected_samples` with your sample names
- [ ] Updated primer sequences (`FWD` and `REV`)  
- [ ] Created `data/metadata_template.csv`
- [ ] Placed FASTQ files in `data/raw/`
- [ ] Successfully ran Step 1 validation
- [ ] Reviewed quality plots from Step 2
- [ ] Adjusted DADA2 parameters if needed

## Common Issues and Solutions

### Issue: "Required FASTQ files are missing"
**Solution:** Check file naming. Files must be exactly `SAMPLENAME_R1.fastq.gz`

### Issue: Low read retention in DADA2
**Solution:** 
1. Check quality plots from Step 2
2. Adjust `truncLen` to less stringent values
3. Increase `maxEE` (e.g., from c(2,2) to c(3,3))

### Issue: "Pipeline functions not found" 
**Solution:** Ensure `Pipeline/Functions/` directory exists with .R files

### Issue: No taxonomic assignments
**Solution:** Run `Rscript scripts/download_databases.R` to get reference databases
