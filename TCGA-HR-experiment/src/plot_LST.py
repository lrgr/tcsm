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
    if row["FANCF"] == 1:
        return "Other HR"
        # return "$\it{ATM}$/$\it{CHEK2}$/$\it{FANCF}$/$\it{FANCM}$"
    if row["FANCM"] == 1:
        return "Other HR"
    if row["ATM"] == 1:
        return "Other HR"
    if row["CHEK2"] == 1:
        return "Other HR"
    return "Wildtype"

if __name__ == '__main__':
    output = []
    for f in snakemake.input:
        output.append(pd.read_csv(f, sep="\t", index_col=0))
    df = pd.concat(output)
    gene = "BRCA1BRCA2RAD51C"
    # print(df["HR Status"])
    heldout = "heldout.ratio"
    df = df.loc[df[gene] == 0]
    print(ranksums(df.loc[df[heldout] > 0]["LST"], df.loc[df[heldout] < 0]["LST"]))
    df["Status of HR Genes"] = df.apply(get_hr_status, axis=1)
    palette ={"Wildtype":"C0","Other HR":"C2"}
    ax = sns.scatterplot(x=df[heldout], y=df["LST"], hue = df.apply(get_hr_status, axis=1), palette=palette)
    ax.set(xlabel="Held-out Log Likelihood Ratio (LLR)", ylabel="LST Count")
    plt.savefig(snakemake.output[0])
