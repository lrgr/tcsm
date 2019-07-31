import pandas as pd
import re
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
#matplotlib.rc('text', usetex = True)

def plot_likelihoods(df, ax):
    for c in df.columns:
        ax.plot(df[c], label=c)
    ax.set_xlabel("Number $\it{K}$ Signatures Used", fontsize=14)
    ax.margins(y=.5)
    ax.legend()
    ax.set_ylabel("Held-out Log Likelihood", fontsize=14)
    return ax


def combine_likelihood_files(files):
    output = []
    for f in files:
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
    return df

if __name__ == '__main__':
    df = combine_likelihood_files(snakemake.input)
    f, ax = plt.subplots()
    ax = plot_likelihoods(df, ax)
    plt.savefig(snakemake.output[0])
