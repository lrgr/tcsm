import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams.update({'font.size': 12})
from textwrap import fill


def get_hr_status(row):
    if row["BRCA1"] == 1:
        return "$\it{BRCA1}$/$\it{BRCA2}$/$\it{RAD51C}$"
    if row["BRCA2"] == 1:
        return "$\it{BRCA1}$/$\it{BRCA2}$/$\it{RAD51C}$"
    if row["RAD51C"] == 1:
        return "$\it{BRCA1}$/$\it{BRCA2}$/$\it{RAD51C}$"
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

def plot_LST(df, axe=None):
    gene = "BRCA1BRCA2RAD51C"
    heldout = "heldout.ratio"
    df["Status of HR Genes"] = df.apply(get_hr_status, axis=1)
    df = df.loc[df["Status of HR Genes"] != "$\it{BRCA1}$/$\it{BRCA2}$/$\it{RAD51C}$"]
    palette ={"Wildtype":"C0","Other HR":"C1", '$\\it{BRCA1}$/$\\it{BRCA2}$/$\\it{RAD51C}$': "C2"}
    ax = sns.scatterplot(x=df[heldout], y=df["LST"], ax=axe, hue = df.apply(get_hr_status, axis=1), palette=palette)
    ax.set(xlabel="Held-out Log Likelihood Ratio (LLR)", ylabel="LST Count")
    return ax

if __name__ == '__main__':
    output = []
    for f in snakemake.input:
        output.append(pd.read_csv(f, sep="\t", index_col=0))
    df = pd.concat(output)
    ax = plot_LST(df)
    ax.plot()
    plt.savefig(snakemake.output[0])
