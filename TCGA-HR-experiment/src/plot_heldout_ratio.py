import pandas as pd
from scipy.stats import ranksums
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams.update({'font.size': 12})
from textwrap import fill
from scipy.stats import ranksums

def get_hr_status(row):
    if row["biBRCA1"] == 1:
        return "biBRCA1"
    if row["biBRCA2"] == 1:
        return "biBRCA2"
    if row["esBRCA1"] == 1:
        return "esBRCA1"
    if row["esRAD51C"] == 1:
        return "esRAD51C"
    if row["esFANCF"] == 1:
        return "esFANCF"
    if row["biFANCM"] == 1:
        return "biFANCM"
    if row["biATM"] == 1:
        return "biATM"
    if row["biCHEK2"] == 1:
        return "biCHEK2"
    return "wildtype"


if __name__ == '__main__':
    output = []
    for f in snakemake.input:
        output.append(pd.read_csv(f, sep="\t", index_col=0))
    df = pd.concat(output)
    genes = list(df.columns)
    gene = "BRCA12"
    heldout = "heldout.ratio"
    df[gene] = df[gene].map({1: "Biallelic Inactivation", 0: "Wildtype"})
    hr_genes = ["biBRCA1", "biBRCA2", "biATM", "biCHEK2", "biFANCM", "esBRCA1", "esRAD51C", "esFANCF"]
    df["HR Status"] = df.apply(get_hr_status, axis=1)
    ax = sns.swarmplot(x=gene, y=heldout, hue="HR Status", data=df)
    # plt.legend(loc='lower center')
    plt.legend(ncol=2)
    ax = sns.boxplot(x=gene, y=heldout, data=df, showcaps=False, boxprops={'facecolor':'None'}, showfliers=False, whiskerprops={'linewidth':0})
    ax.set_xlabel('Held-out Samples')

    plt.savefig(snakemake.output[0])
