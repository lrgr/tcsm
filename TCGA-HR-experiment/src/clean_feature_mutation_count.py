import pandas as pd

if __name__ == '__main__':
    mc_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    bi_inactivation_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    silencing_df = pd.read_csv(snakemake.input[2], sep="\t", index_col=0)
    lst_df = pd.read_csv(snakemake.input[3], sep="\t", header=None, names=["Cancer Type", "Sample", "LST"])
    lst_df = lst_df.drop(columns="Cancer Type")
    lst_df = lst_df.set_index("Sample")
    mc_df.index = [i[0:12] for i in mc_df.index]
    samples = list(set(mc_df.index).intersection(bi_inactivation_df.index).intersection(lst_df.index).intersection(silencing_df.index))
    mc_df = mc_df.loc[samples]
    mc_df.to_csv(snakemake.output[0], sep="\t")
    bi_inactivation_df = bi_inactivation_df.loc[samples]
    lst_df = lst_df.loc[samples]
    silencing_df = silencing_df.loc[samples]
    # now need to combine silencing_df with bi_inactivation_df
    genes = list(set(bi_inactivation_df.columns).intersection(silencing_df.columns))
    silencing_df = silencing_df[genes]
    bi_inactivation_df = bi_inactivation_df[genes]
    feature_df = bi_inactivation_df + silencing_df
    feature_df["BRCA12"] = feature_df["BRCA1"] + feature_df["BRCA2"]
    feature_df = feature_df.merge(lst_df, left_index=True, right_index=True)
    feature_df.to_csv(snakemake.output[1], sep="\t")
