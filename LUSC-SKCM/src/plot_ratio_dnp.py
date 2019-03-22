import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

if __name__ == '__main__':
    df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    heldout = "heldout.ratio"
    feature = "CancerType"
    DNP_cutoff = 15
    # first plot the points with DNP count <= DNP_cutoff
    ax = sns.swarmplot(x=feature, y=heldout, data=df.loc[df["DNPs"]<=DNP_cutoff], marker="o", label="a")
    # now plot the points with DNP count > DNP_cutoff
    ax = sns.swarmplot(x=feature, y=heldout, data=df.loc[df["DNPs"]>DNP_cutoff], marker="v", label="b")
    handles, labels = ax.get_legend_handles_labels()
    labels = ["LUSC <= {} CC>TT".format(DNP_cutoff), "SKCM <= {} CC>TT".format(DNP_cutoff),
              "LUSC > {} CC>TT".format(DNP_cutoff), "SKCM > {} CC>TT".format(DNP_cutoff)]
    ax.legend([handle for i,handle in enumerate(handles)], [label for i,label in enumerate(labels)])
    ax = sns.boxplot(x=feature, y=heldout, data=df, showcaps=False, boxprops={'facecolor':'None'}, showfliers=False, whiskerprops={'linewidth':0})
    ax.set(xlabel="Cancer Type", ylabel="Held-out Log Likelihood Ratio (LLR)")
    plt.savefig(snakemake.output[0])
