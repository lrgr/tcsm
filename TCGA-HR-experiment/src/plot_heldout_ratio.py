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
    if row["FANCM"] == 1:
        return "Other HR"
    if row["ATM"] == 1:
        return "Other HR"
    if row["CHEK2"] == 1:
        return "Other HR"
    return "Wildtype"

def combine_ratio_files(files):
    output = []
    for f in files:
        output.append(pd.read_csv(f, sep="\t", index_col=0))
    return pd.concat(output)

def plot_heldout_ratio(df, axe=None):
    gene = "BRCA1BRCA2RAD51C"
    heldout = "heldout.ratio"
    df[gene] = df[gene].map({1: "Biallelic HR Covariate Inactivation", 0: "Wildtype"})
    df["HR Status"] = df.apply(get_hr_status, axis=1)
    palette ={"Wildtype":"C0","Other HR":"C1", '$\\it{BRCA1}$/$\\it{BRCA2}$/$\\it{RAD51C}$': "C2"}
    ax = sns.swarmplot(x="HR Status", y=heldout, data=df, ax=axe, palette=palette)
    ax = sns.boxplot(x="HR Status", y=heldout, data=df, ax=axe, showcaps=False, boxprops={'facecolor':'None'}, showfliers=False, whiskerprops={'linewidth':0})
    ax.set(xlabel="", ylabel="Held-out Log Likelihood Ratio (LLR)")

if __name__ == '__main__':
    df = combine_ratio_files(snakemake.input)
    ax = plot_heldout_ratio(df)
    ax.plot()
    plt.savefig(snakemake.output[0])
