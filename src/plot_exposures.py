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
    exposure_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    feature_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    df = exposure_df.merge(feature_df, left_index=True, right_index=True)
    topic = snakemake.wildcards["topic"]
    gene=snakemake.wildcards["category"]
    ax = sns.swarmplot(x=gene, y=topic, data=df)
    # plt.legend()
    plt.legend(ncol=1, loc='upper center')
    ax = sns.boxplot(x=gene, y=topic, data=df, showcaps=False, boxprops={'facecolor':'None'}, showfliers=False, whiskerprops={'linewidth':0})
    #ax.set_xlabel('Held-out Samples')
    # plt.show()
    plt.savefig(snakemake.output[0])
