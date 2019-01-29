import pandas as pd
import re
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
matplotlib.rc('text', usetex = True)


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
    df.index = [m.replace("DNPs", "CC$>$TT Count") for m in df.index]
    df.index = [m.replace("Smoker", "Smoking Status") for m in df.index]
    df.index = [m.replace("CancerType", "Cancer Type") for m in df.index]
    df.index = [m.replace("BRCA1BRCA2RAD51C", "$\it{BRCA1}$/$\it{BRCA2}$/${RAD51C}$") for m in df.index]
    df.index = [m.replace("HRINACTIVATIONS", "$\it{ATM}$/$\it{BRCA1}$/$\it{BRCA2}$/$\it{CHEK2}$/$\it{FANCF}$/$\it{FANCM}$/${RAD51C}$") for m in df.index]
    models = df.index.unique()
    df = df.transpose()
    for c in df.columns:
        plt.plot(df[c], label=c)
    plt.xlabel("Number $\it{K}$ Signatures Used")
    plt.legend()
    plt.ylabel("Held-out Log Likelihood")
    # plt.show()
    plt.savefig(snakemake.output[0])
