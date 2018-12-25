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
    heldout = "heldout.ratio"
    print(ranksums(df.loc[df[heldout] > 0]["LST"], df.loc[df[heldout] < 0]["LST"]))
    sns.scatterplot(df[heldout], df["LST"])
    plt.savefig(snakemake.output[0])
