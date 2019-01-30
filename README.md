# Tumor Covariate Signature Model (TCSM)

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

To avoid compatibility issues, we recommend installing the R package `SomaticSignatures`, used for the NMF comparison, in a separate environment.

## Reproducing Figures

### Homologous recombination repair (HR) deficiency in breast cancer

The main code for reproducing these experiments are in the `TCGA-HR-experiment` directory. From the top directory, use the following command to move into the `TCGA-HR-experiment` directory.

```bash
cd TCGA-HR-experiment
```
#### Reproducing Figure 2
```bash
snakemake figure_2
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

#### Reproducing Figure 3

```bash
snakemake figure_3
```

#### Supplemental Figures S4-5

```bash
snakemake supplemental_figures
```
