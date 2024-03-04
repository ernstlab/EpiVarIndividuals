# Epigenetic Variation across Individuals Toy Example

This tutorial demonstrates how to train a 50-state stacked ChromHMM model on a subset of LCL histone mark data (H3K27AC, H3K4ME1, and H2K4ME3) from chromosome 22 across multiple individuals. All steps for running the preprocessing and model training are outlined in the Snakefile. Snakefile can be updated to run on all chromosomes and models with different numbers of states.

## Step 1: Set up environment
### Dependencies
- Snakemake
- Java 1.7 or later
- [ChromHMM](https://compbio.mit.edu/ChromHMM/#:~:text=ChromHMM%20is%20software%20for%20learning,and%20spatial%20patterns%20of%20marks.)
- R
  - Biobase
  - Foreign
  - Mass
  - data.table
### Example conda environment (ChromHMM downloaded separately)
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

## Step 2: Run Snakefile
```
snakemake --cores 1
```
Running the Snakefile will preprocess the signal, binarize the signal for input into ChromHMM, and learn a stacked ChromHMM model.

### Preprocess Signal
- Input: raw histone mark signal and sample covariates
- Output: corrected histone mark signal
- Action: run BA9_reg_cov2_quasipoisson.R
### Binarize Signal
- Input: corrected histone mark signal
- Output: binarized signal
- Action: run ChromHMM BinarizeSignal (see ChromHMM documentation)
### Learn model
- Input: binarized signal
- Output: emission parameters from HMM model and genome segmentation into HMM states
