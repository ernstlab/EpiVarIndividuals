# Genome-wide identification and analysis of recurring patterns of epigenetic variation across individuals

The "stacked" ChromHMM model learns recurring patterns of epigenetic variation across individuals and histone marks throughout the genome. This model was trained on LCL histone mark data published <a href="http://mitra.stanford.edu/kundaje/portal/chromovar3d/">here</a>  and hosted <a href="https://public.hoffman2.idre\linebreak.ucla.edu/ernst/TMZQE/">here</a> for reproducibility. This repository contains code for preprocessing epigenetic data and a toy example demonstrating how to train a stacked ChromHMM model on a subset of the LCL data.

The USCS Genome Browser custom track hub for the LCL genome segmentation and corresponding histone mark data (from the final trained model) can be found <a href="https://public.hoffman2.idre.ucla.edu/ernst/4WSCV/">here</a>.

## Usage of regress_covariates_signal.R for data preprocessing
Covariate correction reduces the number of chromatin states corresponding to sample-level confounders. This process should be performed within each chromosome and mark.

To perform covariate correction on the epigenetic data across individuals, the regress_covariates_signal.R script can be used. This script has 5 positional arguments.
1. A uncorrected signal input file formatted for ChromHMM BinarizeBed 2 A table of covariates with covariates in the columns and samples in the rows
2. Output file for corrected signal
3. Comma separated string of marks (mark ids should be in the sample ids)
4. chromosome id (used for reformatting for ChromHMM BinarizeBed)

## Epigenetic Variation across Individuals Toy Example

This tutorial demonstrates how to train a 50-state stacked ChromHMM model on a subset of LCL histone mark data (H3K27AC, H3K4ME1, and H2K4ME3) from chromosome 22 across multiple individuals. All steps for running the preprocessing and model training are outlined in the Snakefile. Snakefile can be updated to run on all chromosomes and models with different numbers of states.

### Step 1: Set up environment
#### Dependencies
- Snakemake
- Java 1.7 or later
- [ChromHMM](https://compbio.mit.edu/ChromHMM/#:~:text=ChromHMM%20is%20software%20for%20learning,and%20spatial%20patterns%20of%20marks.)
- R
  - Biobase
  - Foreign
  - Mass
  - data.table
#### Example conda environment (ChromHMM downloaded separately)
```
conda create --name stackedhmm_tutorial
conda activate stackedhmm_tutorial
conda install bioconda::snakemake
conda install bioconda::java-jdk
conda install bionconda::bioconductor-biobase
conda install bioconda::bioconductor-biobase
conda install conda-forge::r-foreign
conda install conda-forge::r-mass
conda install conda-forge::r-data.table
conda deactivate
```

### Step 2: Run Snakefile
```
snakemake --cores 1
```
Running the Snakefile will preprocess the signal, binarize the signal for input into ChromHMM, and learn a stacked ChromHMM model.

#### Preprocess Signal
- Input: raw histone mark signal and sample covariates
- Output: corrected histone mark signal
- Action: run BA9_reg_cov2_quasipoisson.R
#### Binarize Signal
- Input: corrected histone mark signal
- Output: binarized signal
- Action: run ChromHMM BinarizeSignal (see ChromHMM documentation)
#### Learn model
- Input: binarized signal
- Output: emission parameters from HMM model and genome segmentation into HMM states
