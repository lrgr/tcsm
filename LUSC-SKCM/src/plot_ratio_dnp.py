import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

# df is a pandas dataframe that contains information about each sample,
# including the "heldout.ratio", the number of "DNPs" and the cancer type of the sample
def plot_likelihood_ratio(df, axe=None):
    heldout = "heldout.ratio"
    feature = "LUSC"
    DNP_cutoff = 15
    # first plot the points with DNP count <= DNP_cutoff
    ax = sns.swarmplot(x=feature, y=heldout, data=df.loc[df["DNPs"]<=DNP_cutoff], marker="^", label="a", ax=axe)
    # now plot the points with DNP count > DNP_cutoff
    ax = sns.swarmplot(x=feature, y=heldout, data=df.loc[df["DNPs"]>DNP_cutoff], marker="v", label="b", ax=axe)
    handles, _ = ax.get_legend_handles_labels()
    labels = ["SKCM <= {} CC>TT".format(DNP_cutoff), "LUSC <= {} CC>TT".format(DNP_cutoff),
              "SKCM > {} CC>TT".format(DNP_cutoff), "LUSC > {} CC>TT".format(DNP_cutoff)]
    ax.legend(handles, labels, loc='lower center')
    ax = sns.boxplot(x=feature, y=heldout, data=df, showcaps=False, boxprops={'facecolor':'None'}, showfliers=False, whiskerprops={'linewidth':0}, ax=axe)
    ax.set_ylabel("Held-out Log Likelihood Ratio (LLR)", fontsize=14)
    ax.set_xlabel("")
    ax.set_xticklabels(["Melanoma (SKCM)", "Lung (LUSC)"])
    ax.tick_params(axis='x', which='major', labelsize=14)
    return ax

if __name__ == '__main__':
    df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    ax = plot_likelihood_ratio(df)
    ax.plot()
    plt.savefig(snakemake.output[0])
