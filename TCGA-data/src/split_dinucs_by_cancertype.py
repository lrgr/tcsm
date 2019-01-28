import pandas as pd

if __name__ == '__main__':
    mc_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    df = pd.read_csv(snakemake.input[1])
    df.index = df.apply(lambda x: "{}>{}".format(x['Ref'], x['Var']), axis=1)
    df = df.drop(columns=['Ref', 'Var'])
    df.columns = [c.split("::")[1] for c in df.columns]
    df = df.transpose()
    # there are some duplicated samples so just take the mean dinucs
    df = df.groupby(df.index).mean()
    print(len(set(mc_df.index)))
    print(len(set(mc_df.index).intersection(df.index)))
    df = df.reindex(mc_df.index)
    df = df[["CC>TT", "GG>AA"]].sum(axis=1).sort_values()
    df.to_csv(snakemake.output[0], sep="\t")
