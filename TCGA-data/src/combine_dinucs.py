import pandas as pd

if __name__ == '__main__':
    alex_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0, names=["count"])
    alex_df.index = [i[0:12] for i in alex_df.index]
    firehose_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    firehose_df.index = [i[0:12] for i in firehose_df.index]
    firehose_df = firehose_df.loc[firehose_df["Category"].isin(["CC>TT", "GG>AA"])]
    firehose_df = firehose_df.groupby(firehose_df.index).sum()
    print(alex_df.shape)
    print(firehose_df.shape)
    df = pd.concat([alex_df, firehose_df])
    df = df.groupby(df.index).mean()
    df.columns=["DNPs"]
    df.to_csv(snakemake.output[0], sep="\t")
