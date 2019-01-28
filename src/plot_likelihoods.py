import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams['font.family'] = 'Courier New'


if __name__ == '__main__':
    output = []
    for f in snakemake.input:
        # get the covariates used from the filename
        out = re.split("_|\.", f)
        covariates = out[1]
        # stored as one column with one row for each K
        df = pd.read_csv(f, sep="\t", index_col=0).transpose()
        df.index = [covariates]
        output.append(df)
    df = pd.concat(output)
    df.index = [m.replace("NULL", "No Covariates") for m in df.index]
    df.index = [m.replace("BRCA1BRCA2RAD51C", "Biallelic HR Covariate") for m in df.index]
    df.index = [m.replace("HRINACTIVATIONS", "Biallelic Inactivtion of any HR Gene") for m in df.index]
    models = df.index.unique()
    df = df.transpose()
    for c in df.columns:
        plt.plot(df[c], label=c)
    plt.xlabel("Number K Signatures Used")
    plt.legend()
    plt.ylabel("Log Likelihood")
    # plt.show()
    plt.savefig(snakemake.output[0])
