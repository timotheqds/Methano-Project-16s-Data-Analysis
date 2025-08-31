# MethanoPipeline Troubleshooting Guide

## Installation Issues

### R Package Installation Fails

**Error:** Package installation timeouts or compilation errors

**Solutions:**
```bash
# Option 1: Use binary packages (faster)
Rscript -e "install.packages('BiocManager', type='binary')"

# Option 2: Install packages individually
Rscript -e "BiocManager::install('dada2')"
Rscript -e "install.packages('tidyverse')"

# Option 3: Use conda for some packages
conda install -c bioconda r-dada2
```

### Conda Environment Issues

**Error:** `conda env create` fails

**Solutions:**
```bash
# Update conda first
conda update conda

# Create environment with specific channels
conda create -n methanopipeline -c bioconda -c conda-forge r-base=4.3 cutadapt fastqc

# Activate and install remaining packages
conda activate methanopipeline
pip install lulu
```

## Data Validation Issues

### "Required FASTQ files are missing"

**Check these:**
1. File naming convention: `SAMPLENAME_R1.fastq.gz` and `SAMPLENAME_R2.fastq.gz`
2. Files are in correct directory: `data/raw/`
3. Sample names in script match your files exactly

**Debug commands:**
```bash
# List your files
ls data/raw/*.fastq.gz

# Check naming pattern
ls data/raw/ | grep -E "_R[12]\.fastq\.gz$"
```

### "Pipeline functions not found"

**Solutions:**
1. Ensure `Pipeline/Functions/` directory exists
2. Check if .R files are present: `ls Pipeline/Functions/*.R`
3. If missing, you need the original SimpleMetaPipeline functions

## Quality Control Issues

### Low Read Retention After Filtering

**Symptoms:** Very few reads pass DADA2 filtering

**Solutions:**
1. **Check quality plots** from Step 2:
   ```bash
   # Look at quality profiles
   open data/results/quality_profiles_forward.pdf
   open data/results/quality_profiles_reverse.pdf
   ```

2. **Adjust parameters** in Step 1:
   ```r
   # More lenient parameters
   truncLen = c(200, 150)    # Reduce truncation lengths
   maxEE = c(3, 3)          # Allow more expected errors
   ```

3. **Check primer removal:** Ensure cutadapt worked properly

### Poor Quality Profiles

**Symptoms:** Quality drops very early in reads

**Solutions:**
1. **Check for adapter contamination**
2. **Verify primer sequences** are correct
3. **Consider different truncation strategy**

## DADA2 Processing Issues

### "No reads passed the filter"

**Debug steps:**
1. Check if files exist after cutadapt
2. Verify file paths are correct
3. Try more lenient filtering parameters

**Emergency parameters:**
```r
# Very lenient for testing
truncLen = c(150, 120)
maxEE = c(5, 5)
truncQ = 1
```

### Memory Issues

**Error:** R runs out of memory

**Solutions:**
```r
# Reduce memory usage
multithread = FALSE  # Disable multithreading
pool = FALSE         # Disable pooling
```

## Clustering and Taxonomy Issues

### LULU Fails

**Error:** LULU curation fails

**Solutions:**
1. Check if ASV table has sufficient data
2. Verify match list was created properly
3. Try with different LULU parameters:
   ```r
   MatchRate1 <- 90      # More stringent
   MinRelativeCo1 <- 0.5 # Less stringent
   ```

### No Taxonomic Assignments

**Causes:**
1. Reference database not downloaded
2. Database file path incorrect
3. Sequences too short/poor quality

**Solutions:**
```bash
# Download databases
Rscript scripts/download_databases.R

# Check database exists
ls databases/

# Skip taxonomy for testing
# (Step 4 will continue without taxonomy)
```

## File Path and Permission Issues

### "Permission denied" Errors

**Solutions:**
```bash
# Fix file permissions
chmod -R 755 data/
chmod -R 755 scripts/

# Ensure directories exist
mkdir -p data/{raw,processed,results,step_outputs}
mkdir -p logs/steps
```

### Absolute vs Relative Paths

**Issue:** Scripts fail when run from different directories

**Solution:** Always run from project root:
```bash
cd /path/to/Methano-Project-16s-Data-Analysis
Rscript scripts/steps/step1_setup_validation.R
```

## Step-Specific Troubleshooting

### Step 1 Issues
- Verify sample names in script match your files
- Check metadata file format
- Ensure Pipeline/Functions/ exists

### Step 2 Issues  
- Verify cutadapt is installed: `which cutadapt`
- Check primer sequences are correct
- Ensure sufficient disk space for processed files

### Step 3 Issues
- Review quality profiles before running
- Adjust truncLen based on quality plots
- Check for sufficient reads passing filters

### Step 4 Issues
- Verify previous steps completed successfully
- Check if databases are needed for your analysis
- Monitor memory usage during clustering

## Performance Optimization

### Speed Up Processing

1. **Use multithread:** `multithread = TRUE` (default)
2. **Adjust memory:** Increase R memory limit
3. **Use subset for testing:** Test with 2-4 samples first

### Reduce Resource Usage

```r
# Lighter processing
multithread = FALSE
pool = FALSE
ClusteringTechnique = "vsearch"  # Instead of swarm
```

## Getting Help

### Debug Information
Add to any script for detailed error tracking:
```r
options(error = function() {
  traceback()
  q(status = 1)
})
```

### Log Files
Check execution logs:
```bash
# If using step_controller.R
tail -f logs/steps/step1.log
```

### Create Reproducible Example
For help requests, include:
1. Your modified `expected_samples` vector
2. Your primer sequences  
3. First few lines of `ls data/raw/`
4. Full error message
5. Output of `Rscript scripts/verify_installation.R`

## Emergency Reset

If things get completely broken:
```bash
# Clean all outputs
rm -rf data/processed/
rm -rf data/results/
rm -rf data/step_outputs/

# Recreate directories
mkdir -p data/{processed,results,step_outputs}

# Start over from Step 1
Rscript scripts/steps/step1_setup_validation.R
```
