import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

if __name__ == '__main__':
    df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    #dnp_df = pd.concat([pd.read_csv(snakemake.input[1], sep="\t", index_col=0), pd.read_csv(snakemake.input[2], sep="\t", index_col=0)])
    heldout = "heldout.ratio"
    feature = "CancerType"
    # df = dnp_df.merge(ratio_df, left_index=True, right_index=True)
    print(df.shape)
    DNP_cutoff = 15
    # first plot the points with DNP count <= DNP_cutoff
    ax = sns.swarmplot(x=feature, y=heldout, data=df.loc[df["DNPs"]<=DNP_cutoff], marker=".", label="a")#color="b")
    # now plot the points with DNP count > DNP_cutoff
    ax = sns.swarmplot(x=feature, y=heldout, data=df.loc[df["DNPs"]>DNP_cutoff], marker="*", label="b") #color="r")
    handles, labels = ax.get_legend_handles_labels()
    labels = ["LUSC <= {} CC>TT".format(DNP_cutoff), "SKCM <= {} CC>TT".format(DNP_cutoff),
              "LUSC > {} CC>TT".format(DNP_cutoff), "SKCM > {} CC>TT".format(DNP_cutoff)]
    # handles[0]
    ax.legend([handle for i,handle in enumerate(handles)], [label for i,label in enumerate(labels)])
    # print(ax.legend())
    # m = markers.MarkerStyle(marker="*", label="something")
    # ax = sns.swarmplot(x=feature, y=heldout, data=df)
    # plt.legend(loc='lower center')
    # plt.legend(ncol=2)
    ax = sns.boxplot(x=feature, y=heldout, data=df, showcaps=False, boxprops={'facecolor':'None'}, showfliers=False, whiskerprops={'linewidth':0})
    ax.set(xlabel="Cancer Type", ylabel="Held-out Log Likelihood Ratio (LLR)")
    # plt.legend(handles=[m])
    #plt.colorbar()
    # lusc_df = df.loc[df["CancerType"] == "LUSC"]
    # ax = sns.boxplot(x=lusc_df["heldout.ratio"]>0, y=lusc_df["count"])
    # ax.set_title("LUSC samples")
    # ax.set(xlabel="Predicted to be LUSC", ylabel="Number of CC>TT DNPs")
    # plt.ylim(-5, 200)
    # plt.show()
    # skcm_df = df.loc[df["CancerType"] == "SKCM"]
    # ax = sns.boxplot(x=skcm_df["heldout.ratio"]>0, y=skcm_df["count"])
    # ax.set_title("SKCM samples")
    # ax.set(xlabel="Predicted to be LUSC", ylabel="Number of CC>TT DNPs")
    # plt.ylim(-5, 200)
    # plt.show()
    # print(df)
    # ax = sns.scatterplot(x="heldout.ratio", y="count", data=df, hue="CancerType", hue_order=["SKCM", "LUSC"])
    # ax.set(xlabel="Held-out log likelihood ratio", ylabel="Number of CC>TT DNPs")
    # plt.ylim(-5, 200)
    # plt.show()
    plt.savefig(snakemake.output[0])
