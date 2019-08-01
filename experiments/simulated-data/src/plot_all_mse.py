import pandas as pd
import re
from textwrap import fill
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams.update({'font.size': 12})
rcParams['axes.xmargin'] = 0

def combine_mse_files(files):
    output = []
    for i in range(0, len(files)):
        # extract the number of samples from the filename
        model = re.split("_|\.", snakemake.input[i])[1]
        df = pd.read_csv(snakemake.input[i], sep="\t", index_col=0)
        df["model"] = model
        output.append(df)
    df = pd.concat(output)
    df["model"] = df["model"].map({"stm": "TCSM", "ctm": "TCSM (no covariates)", "SS": "Somatic Signatures" })
    df = df.set_index('model', append=True)
    df.columns = [int(col) for col in df.columns]
    return df.reindex(sorted(df.columns), axis=1)

def plot_mse(df, signature, ax):
    model_color_map = {"TCSM": "r", "TCSM (no covariates)": "b", "Somatic Signatures": "g"}
    signature_df = df.loc[signature]
    models = signature_df.index
    for model in models:
        ax.plot(signature_df.loc[model], c=model_color_map[model], marker="o")
    ax.set_title("COSMIC Signature {}".format(signature[-1]))

if __name__ == '__main__':
    df = combine_mse_files(snakemake.input)
    signatures = ["SBS3", "SBS5"]
    fig, axs = plt.subplots(ncols=len(signatures))
    fig.set_figheight(3)
    fig.set_figwidth(7.5)
    for signature, ax in zip(signatures, axs):
        plot_mse(df, signature, ax)
    axs[0].set_ylabel('Mean Squared Error')
    plt.savefig(snakemake.output[0])
