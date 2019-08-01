import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams.update({'font.size': 12})
rcParams['axes.xmargin'] = 0


if __name__ == '__main__':
    output = []
    for i in range(0, len(snakemake.input)):
        # extract the number of samples from the filename
        model = re.split("_|\.", snakemake.input[i])[1]
        df = pd.read_csv(snakemake.input[i], sep="\t", index_col=0)
        df["model"] = model
        output.append(df)
    df = pd.concat(output)
    df["model"] = df["model"].map({"stm": "TCSM", "ctm": "TCSM (no covariates)", "SS": "Somatic Signatures" })
    df = df.set_index('model', append=True)
    df.columns = [int(col) for col in df.columns]
    df = df.reindex(sorted(df.columns), axis=1)
    # index is the number of samples
    # columns are multi-index, the signature and models (ex. stm and ctm)
    # we want to create a subplot for each signature
    signatures = ["SBS3", "SBS5"]
    model_color_map = {"TCSM": "r", "TCSM (no covariates)": "b", "Somatic Signatures": "g"}
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
    ax[0].set_ylabel('Cosine Similarity')
    # handles, labels = row.get_legend_handles_labels()
    # fig.legend(handles, labels, loc='upper right')
    # plt.show()
    plt.savefig(snakemake.output[0])
