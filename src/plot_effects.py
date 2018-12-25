import pandas as pd
import seaborn as sns


if __name__ == '__main__':
    df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    features = list(df.columns)
    features.remove("(Intercept)")
    for f in features:
        df[f] = df["(Intercept)"] + df[f]
    ax = sns.heatmap(df, annot=True) #, annot_kws={"size": 10})
    fig = ax.get_figure()
    fig.savefig(snakemake.output[0])
