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
    if row["BRCA1"] == 1:
        return "BRCA1"
    if row["BRCA2"] == 1:
        return "BRCA2"
    if row["RAD51C"] == 1:
        return "RAD51C"
    if row["FANCF"] == 1:
        return "FANCF"
    if row["FANCM"] == 1:
        return "FANCM"
    if row["ATM"] == 1:
        return "ATM"
    if row["CHEK2"] == 1:
        return "CHEK2"
    return "wildtype"


if __name__ == '__main__':
    output = []
    for f in snakemake.input:
        output.append(pd.read_csv(f, sep="\t", index_col=0))
    df = pd.concat(output)
    gene = "BRCA1BRCA2RAD51C"
    heldout = "heldout.ratio"
    df[gene] = df[gene].map({1: "Biallelic HR Covariate Inactivation", 0: "Wildtype"})
    df["HR Status"] = df.apply(get_hr_status, axis=1)
    print(ranksums(df.loc[df[gene] == "Biallelic HR Covariate Inactivation"][heldout], df.loc[df[gene] == "Wildtype"][heldout]))
    ax = sns.swarmplot(x=gene, y=heldout, hue="HR Status", data=df)
    # plt.legend(loc='lower center')
    plt.legend(ncol=2)
    ax = sns.boxplot(x=gene, y=heldout, data=df, showcaps=False, boxprops={'facecolor':'None'}, showfliers=False, whiskerprops={'linewidth':0})
    ax.set(xlabel="", ylabel="Held-out log likelihood ratio")

    plt.savefig(snakemake.output[0])
