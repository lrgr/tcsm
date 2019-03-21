import pandas as pd
import re
from textwrap import fill
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams.update({'font.size': 12})
rcParams['axes.xmargin'] = 0


if __name__ == '__main__':
    output = []
    for f in snakemake.input:
        # extract the number of samples from the filename
        in_file = re.split("_|\.", f)
        covariates = in_file[2]
        model = in_file[1]
        df = pd.read_csv(f, sep="\t", index_col=0)
        df["model"] = "{}_{}".format(model, covariates)
        output.append(df)
    df = pd.concat(output)
    df["model"] = df["model"].map({"stm_feature1": "TCSM", "stm_NULL": "TCSM (no covariates)", "SS_NULL": "Somatic Signatures" })
    df = df.set_index('model', append=True)
    df.columns = [int(col) for col in df.columns]
    df = df.reindex(sorted(df.columns), axis=1)
    # index is the number of samples
    # columns are multi-index, the signature and models (ex. stm and ctm)
    # we want to create a subplot for each signature
    signatures = ["SBS3", "SBS5"]
    model_color_map = {"TCSM": "r", "TCSM (no covariates)": "b", "Somatic Signatures": "g"}
    # print(df)
    fig, ax = plt.subplots(ncols=len(signatures))
    fig.set_figheight(3)
    fig.set_figwidth(7.5)
    for signature, row in zip(signatures, ax):
        signature_df = df.loc[signature]
        models = signature_df.index
        for model in models:
            row.plot(signature_df.loc[model], c=model_color_map[model], marker="o")
        row.set_title("COSMIC Signature {}".format(signature[-1]))
        # row.set_xlabel("Number of Samples")
    ax[0].set_ylabel('Mean Squared Error')
    handles, labels = row.get_legend_handles_labels()
    #labels = [fill(l, 8, break_long_words=False) for l in labels]
    #row.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.04,1.1))
    # fig.legend(handles, labels, loc='upper right')

    # fig.text(0.0, 0.5, 'mean MSE', va='center', rotation='vertical')
    # plt.xlabel("number of samples")
    # plt.show()
    plt.savefig(snakemake.output[0])
