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
    gene = "SKCM"
    heldout = "heldout.ratio"
    df[gene] = df[gene].map({1: "Melanoma", 0: "LUSC"})
    print(df.loc[df[gene] == "LUSC"].sort_values(by=heldout))
    #print(ranksums(df.loc[df[gene] == "Biallelic HR Covariate Inactivation"][heldout], df.loc[df[gene] == "Wildtype"][heldout]))
    ax = sns.swarmplot(x=gene, y=heldout, data=df)
    # plt.legend(loc='lower center')
    plt.legend(ncol=2)
    ax = sns.boxplot(x=gene, y=heldout, data=df, showcaps=False, boxprops={'facecolor':'None'}, showfliers=False, whiskerprops={'linewidth':0})
    ax.set(xlabel="", ylabel="Held-out log likelihood ratio")

    plt.savefig(snakemake.output[0])
