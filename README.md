# Tumor Covariate Signature Model (TCSM)
<img src='https://travis-ci.org/lrgr/tcsm.svg?branch=master'>

This repository contains the code for reproducing the results on real data from the paper Modeling Clinical and Molecular Covariates of Mutational Process Activity in Cancer, which implements the Tumor Covariate Signature Model (TCSM). The code for reproducing the simulated data results from the same paper are available upon request. TCSM requires Python 3. We recommend using Anaconda to install dependencies.

Once you have Anaconda installed, you can create an environment with the dependencies and activate it with the following commands:
```bash
conda env create -f envs/environment.yml
conda activate tcsm
```

TCSM requires two submodules, mutation-signature-viz and signature-estimation-py, which can be installed by the following commands
```
git submodule init
git submodule update
```

## References

Welles Robinson, Roded Sharan, Max Leiserson. (2019). Modeling clinical and molecular covariates of mutational process activity in cancer. _Bioinformatics_. [paper link](https://doi.org/10.1093/bioinformatics/btz340)

## Demo

We have provided a demo (`demo`) to help users run TCSM on their own datasets. The demo shows how to use TCSM on real data with and without using covariates. First, the demo estimates the number of signatures in the dataset by plotting the heldout log-likelihood across a range of K. Next, we show how to estimate signatures and exposures using TCSM with and without covariates. Finally, for an advanced use case, we should how to estimate exposures with TCSM when the covariate value is unknown or hidden for a subset of the samples.

### Key Commands

The key command is `src/run_stm.R`, which is used to estimate signatures and exposures on a dataset. It requires a mutation count matrix and the number of signatures to use. It outputs signatures, exposures and correlations between the signatures (sigma).

If covariates are used, it expects covariate values for all samples and the names of covariates to use (separated by +). In addition to the above, it outputs an estimation of the effect of each covariate on each signature (effect table) as well as the parameters to the logistic normal (gamma), which is another (difficult to directly interpret) measure of the effect of each covariate on each signature.

## Reproducing Key Results

### Homologous recombination repair (HR) deficiency in breast cancer

We use Travis CI to regenerate Figure 3 every time the code base is updated.

<!-- <img src='http://tcsm.lrgr.io/fig2.png'>

**Figure 2**: Benchmark of TCSM with (red) and without (blue) covariates and NMF-based SomaticSignatures (green) on synthetic data. (A) Cosine similarity of inferred signatures (β) to hidden Signatures 3 and 5 using the true K = 4 averaged across 50 datasets, varying the number of samples. (B) Mean-squared error of the inferred exposures (θ) for the same datasets as in (A) -->

<img src='http://tcsm.lrgr.io/fig3.png'>

**Figure 3**: (A) Comparison of the log-likelihood of held-out samples across K = 2–10 between TCSM with the biallelic HR covariate (inactivations of BRCA1, BRCA2 or RAD51C) and TCSM without covariates. (B) The log-likelihood ratio (LLR) of samples with the biallelic HR covariate hidden where LLR>0 indicates the mutations of a sample are more likely under the biallelic HR covariate inactivation model. (C) After excluding tumors with known biallelic inactivations in BRCA1, BRCA2 or RAD51C, the plot of a tumor’s LLR against its LST count



To reproduce this figure and other experiments related to HR-deficiency, first go to the `TCGA-HR-experiment` directory and then use snakemake.

```bash
cd TCGA-HR-experiment
snakemake figure_3
```

### Simultaneously learning signatures in melanomas and lung cancers
The main code for reproducing these experiments are in the `LUSC-SKCM` directory. From the top directory, use the following command to move into the `LUSC-SKCM` directory.

```bash
cd LUSC-SKCM
```

To reproduce Figure 4, use the following command.

```bash
snakemake figure_4
```
