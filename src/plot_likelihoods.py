import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams['font.family'] = 'Courier New'


if __name__ == '__main__':
    output = []
    for i in range(0, len(snakemake.input)):
        out = re.split("_|\.", snakemake.input[i])
        model = out[1]
        print(model)
        df = pd.read_csv(snakemake.input[i], sep="\t", index_col=0).transpose()
        df.index = [model]
        output.append(df)
    df = pd.concat(output)
    #n_samples_groups = df["n_samples"].unique()
    print(df.index)
    #df.index = df.index.map({"stm": "TCSM", "ctm": "TCSM (no covariates)" })



    plt.plot(df.loc["stm"].transpose(), color="r", marker="o")
    plt.plot(df.loc["ctm"].transpose(), color="b", marker="o")
    plt.xlabel("Number K Signatures Used")
    plt.legend()
    plt.ylabel("Log Likelihood")
    # plt.show()
    # fig is a Figure object
    # ax is an array of Axes objects
    # fig, ax = plt.subplots(ncols=len(n_samples_groups), nrows=1)
    # fig.set_figheight(3)
    # fig.set_figwidth(15)
    # for n_sample_group, row in zip(n_samples_groups, ax):
    #     sample_df = df.loc[df["n_samples"] == n_sample_group]
    #     sample_df = sample_df.drop(columns="n_samples")
    #     row.plot(sample_df.loc["TCSM"].transpose(), color="r", marker="o")
    #     row.plot(sample_df.loc["TCSM (no covariates)"].transpose(), color="b", marker="o")
    #     # row.legend()
    #     row.set_title("{} samples".format(n_sample_group))
    #
    #     row.set_xlabel("Number K Signatures Used")
    # handles, labels = row.get_legend_handles_labels()
    # fig.legend(handles, labels, loc='upper right')
    # ax[0].set_ylabel("Held-out Log Likelihood")
    #     # plt.xticks([sample_df.index])
    #     # row.set_title("{} samples".format(n_sample_group), x=-0.1, y=0.5)
    # # row.xlabel("number of signatures used")
    # # fig.text(0.0, 0.5, 'held-out log likelihood', va='center', rotation='vertical')
    # #
    # # fig.tight_layout()
    # # plt.tight_layout(pad=1.2)
    # # plt.xlabel("number of signatures used")
    # # plt.legend()
    # plt.show()
    plt.savefig(snakemake.output[0])
