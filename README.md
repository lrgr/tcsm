# Tumor Covariate Signature Model (TCSM)
<img src='https://travis-ci.org/lrgr/tcsm.svg?branch=master'>

This repository contains the code for reproducing the results on real data from the paper Modeling Clinical and Molecular Covariates of Mutational Process Activity in Cancer, which implements the Tumor Covariate Signature Model (TCSM).

## References

Welles Robinson, Roded Sharan, Max Leiserson. (2019). Modeling clinical and molecular covariates of mutational process activity in cancer. _Bioinformatics_. [paper link](https://doi.org/10.1093/bioinformatics/btz340)

## Basic Usage

TCSM requires Python 3. We recommend using Anaconda to install dependencies.

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

### Key Commands

To run TCSM without covariates, use the following command
```
./src/run_tcsm.R <mutation-count-input-file> <num-signatures>
```

To run TCSM with covariates, use the following command
```
./src/run_tcsm.R <mutation-count-input-file> <num-signatures> -c=<covariate-file-input> --covariates <covariate-names(separated by +)>
```
### Demo

We have provided a demo (`demo`) to help users run TCSM on their own datasets. The demo shows how to use TCSM on real data with and without using covariates. First, the demo estimates the number of signatures in the dataset by plotting the heldout log-likelihood across a range of K. Next, we show how to estimate signatures and exposures using TCSM with and without covariates. Finally, for an advanced use case, we should how to estimate exposures with TCSM when the covariate value is unknown or hidden for a subset of the samples.

## Reproducing Key Results

We use Travis CI to regenerate the key figures of the paper (Figure 3 and 4) when the master branch is updated.

### Homologous recombination repair (HR) deficiency in breast cancer

<img src='http://lrgr.io/tcsm/figure3.png'>

**Figure 3**: (A) Comparison of the log-likelihood of held-out samples across K = 2–10 between TCSM with the biallelic HR covariate (inactivations of BRCA1, BRCA2 or RAD51C) and TCSM without covariates. (B) The log-likelihood ratio (LLR) of samples with the biallelic HR covariate hidden where LLR>0 indicates the mutations of a sample are more likely under the biallelic HR covariate inactivation model. (C) After excluding tumors with known biallelic inactivations in BRCA1, BRCA2 or RAD51C, the plot of a tumor’s LLR against its LST count


### Simultaneously learning signatures in melanomas and lung cancers

<img src='http://lrgr.io/tcsm/figure4.png'>

**Figure 4**: (A) The heldout log-likelihood plot used for model selection to obtain K = 4. (B) The log-likelihood ratio (LLR) of the cancer type covariate for tumors where LLR <0 means the mutations of the tumor are more likely under LUSC and LLR >0 means the mutations of the tumor are more likely under SKCM
