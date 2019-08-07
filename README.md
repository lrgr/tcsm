# Tumor Covariate Signature Model (TCSM)
<img src='https://travis-ci.org/lrgr/tcsm.svg?branch=master'>

This repository contains the code for reproducing the results on real data from the paper Modeling Clinical and Molecular Covariates of Mutational Process Activity in Cancer, which implements the Tumor Covariate Signature Model (TCSM). The code for reproducing the simulated data results from the same paper are available upon request. TCSM requires Python 3. We recommend using Anaconda to install dependencies.

Once you have Anaconda installed, you can create an environment with the dependencies and activate it with the following commands:
```bash
conda env create -f environment.yml
conda activate tcsm
```

TCSM requires two submodules, mutation-signature-viz and signature-estimation-py, which can be installed by the following commands
```
git submodule init
git submodule update
```

## References

Welles Robinson, Roded Sharan, Max Leiserson. (2019). Modeling clinical and molecular covariates of mutational process activity in cancer. _Bioinformatics_. [paper link](https://doi.org/10.1093/bioinformatics/btz340)

## Figures

<!-- <img src='http://tcsm.lrgr.io/fig2.png'>

**Figure 2**: Benchmark of TCSM with (red) and without (blue) covariates and NMF-based SomaticSignatures (green) on synthetic data. (A) Cosine similarity of inferred signatures (β) to hidden Signatures 3 and 5 using the true K = 4 averaged across 50 datasets, varying the number of samples. (B) Mean-squared error of the inferred exposures (θ) for the same datasets as in (A) -->

<img src='http://lrgr.io/tcsm/figure3.png'>

**Figure 3**: (A) Comparison of the log-likelihood of held-out samples across K = 2–10 between TCSM with the biallelic HR covariate (inactivations of BRCA1, BRCA2 or RAD51C) and TCSM without covariates. (B) The log-likelihood ratio (LLR) of samples with the biallelic HR covariate hidden where LLR>0 indicates the mutations of a sample are more likely under the biallelic HR covariate inactivation model. (C) After excluding tumors with known biallelic inactivations in BRCA1, BRCA2 or RAD51C, the plot of a tumor’s LLR against its LST count

### Homologous recombination repair (HR) deficiency in breast cancer

The main code for reproducing these experiments are in the `TCGA-HR-experiment` directory. From the top directory, use the following command to move into the `TCGA-HR-experiment` directory.

```bash
cd TCGA-HR-experiment
```
#### Reproducing Figure 3
```bash
snakemake figure_3
```

#### Reproducing Biallelic Inactivation Prediction

To calculate the AUPRC using the exposures estimated by TCSM with the biallelic covariate and TCSM without covariates
```bash
eval_HR_prediction_stm
```

To calculate the AUPRC using the exposures estimated by SomaticSignatures (should be in separate environment with SomaticSignatures installed)
```bash
eval_HR_prediction_SS
```

#### Reproducing Supplemental Figures S1-3

To reproduce figures S1 and S2
```bash
snakemake supplemental_figures_S1_S2
```

To reproduce figure S3
```bash
snakemake supplemental_figures_S3
```

### Simultaneously learning signatures in melanomas and lung cancers
The main code for reproducing these experiments are in the `LUSC-SKCM` directory. From the top directory, use the following command to move into the `LUSC-SKCM` directory.

```bash
cd LUSC-SKCM
```

#### Reproducing Figure 4

```bash
snakemake figure_4
```

#### Supplemental Figures S4-5

```bash
snakemake supplemental_figures
```
