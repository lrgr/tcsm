import pandas as pd
from scipy.stats import ranksums
from statistics import mean

if __name__ == '__main__':
    df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    dnp_df = pd.concat([pd.read_csv(snakemake.input[1], sep="\t", index_col=0), pd.read_csv(snakemake.input[2], sep="\t", index_col=0)])
    df = dnp_df.merge(df, left_index=True, right_index=True)
    df = df.loc[df["CancerType"] == "SKCM"]
    print(mean(df.loc[df["heldout.ratio"] > 0]["nMuts"]))
    print(mean(df.loc[df["heldout.ratio"] < 0]["nMuts"]))
    print(ranksums(df.loc[df["heldout.ratio"] > 0]["nMuts"], df.loc[df["heldout.ratio"] < 0]["nMuts"]))
    print(mean(df.loc[df["heldout.ratio"] > 0]["count"]))
    print(mean(df.loc[df["heldout.ratio"] < 0]["count"]))
    print(ranksums(df.loc[df["heldout.ratio"] > 0]["count"], df.loc[df["heldout.ratio"] < 0]["count"]))
