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


if __name__ == '__main__':
    output = []
    for f in snakemake.input:
        output.append(pd.read_csv(f, sep="\t", index_col=0))
    df = pd.concat(output)
    genes = list(df.columns)
    # genes.remove("heldout.ratio")
    # genes.remove("BRCA1")
    # genes.remove("BRCA2")
    # genes.remove("BRCA12")
    # print(df[genes].sum(axis=1))
    #df["HRd"] = (df[genes].sum(axis=1) > 0).astype(dtype="int32")
    #print(df.sum(axis=0).sort_values())
    gene = "BRCA12"
    heldout = "heldout.ratio"
    # df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    #heldout = "Log-Likelihood Ratio of $\itBRCA$12 Inactivated versus Wildtype Models"
    #heldout = fill(heldout, 40, break_long_words=False)
    #df = df.rename(columns={"BRCA12": "BRCA Status", "heldout.ratio": heldout})
    df[gene] = df[gene].map({1: "Biallelic Inactivation", 0: "Wildtype"})
    # df = df.loc[df[gene] == "Wildtype"]
    print(ranksums(df.loc[df[heldout] > 0]["LST"], df.loc[df[heldout] < 0]["LST"]))
    #sns.scatterplot(df[heldout], df["LST"])
    #plt.show()
    # df[heldout]
    # print(df.groupby(by="true.brcaness").mean())
    # sig = ranksums(df.loc[df["BRCA Status"] == "BRCA-deficient"]["Held-out Log-Likelihood Ratio"], df.loc[df["BRCA Status"] == "BRCA-proficient"]["Held-out Log-Likelihood Ratio"])
    # print(sig)
    ax = sns.swarmplot(x=gene, y=heldout, data=df)
    ax = sns.boxplot(x=gene, y=heldout, data=df, showcaps=False, boxprops={'facecolor':'None'}, showfliers=False, whiskerprops={'linewidth':0})
    ax.set_xlabel('Held-out Samples')
    # plt.show()
    plt.savefig(snakemake.output[0])
