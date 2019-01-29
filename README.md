# Tumor Covariate Signature Model (TCSM)

TCSM requires Python 3. We recommend using Anaconda to install dependencies.

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

#### Homologous Recombination-related Figures

```bash
cd TCGA-HR-experiment
snakemake figure_2
```

#### SKCM vs. LUSC comparison Figures

```bash
cd LUSC-SKCM
snakemake figure_3
```
