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
    exposure_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    feature_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    df = exposure_df.merge(feature_df, left_index=True, right_index=True)
    topic = snakemake.wildcards["topic"]
    df["HR Status"] = df.apply(get_hr_status, axis=1)
    gene="BRCA1"
    ax = sns.swarmplot(x=gene, y=topic, hue="HR Status", data=df)
    # plt.legend()
    plt.legend(ncol=1, loc='upper center')
    ax = sns.boxplot(x=gene, y=topic, data=df, showcaps=False, boxprops={'facecolor':'None'}, showfliers=False, whiskerprops={'linewidth':0})
    #ax.set_xlabel('Held-out Samples')
    # plt.show()
    plt.savefig(snakemake.output[0])
