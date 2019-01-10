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
    gene = "BRCA1BRCA2RAD51C"
    # print(df["HR Status"])
    heldout = "heldout.ratio"
    df = df.loc[df[gene] == 0]
    print(ranksums(df.loc[df[heldout] > 0]["LST"], df.loc[df[heldout] < 0]["LST"]))
    ax = sns.scatterplot(x=heldout, y="LST", data=df)
    ax.set(xlabel="Held-out log likelihood ratio", ylabel="Number of LST")
    plt.savefig(snakemake.output[0])
