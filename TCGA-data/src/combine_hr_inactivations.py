import pandas as pd

if __name__ == '__main__':
    bi_inactivation_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    silencing_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    samples = list(set(bi_inactivation_df.index).intersection(silencing_df.index))
    bi_inactivation_df = bi_inactivation_df.loc[samples]
    silencing_df = silencing_df.loc[samples]
    # now need to combine silencing_df with bi_inactivation_df
    genes = list(set(bi_inactivation_df.columns).intersection(silencing_df.columns))
    silencing_df = silencing_df[genes]
    bi_inactivation_df = bi_inactivation_df[genes]
    feature_df = bi_inactivation_df + silencing_df
    feature_df["BRCA1&BRCA2"] = feature_df["BRCA1"] + feature_df["BRCA2"]
    feature_df["BRCA12"] = feature_df["BRCA1"] + feature_df["BRCA2"]
    # feature_df["BRCA1&BRCA2&RAD51C"] = feature_df["BRCA1&BRCA2"] + feature_df["RAD51C"]
    # feature_df["BRCA1&BRCA2&RAD51C&ATM"] = feature_df["BRCA1&BRCA2&RAD51C"] + feature_df["ATM"]
    # feature_df["BRCA1&BRCA2&RAD51C&FANCM"] = feature_df["BRCA1&BRCA2&RAD51C"] + feature_df["FANCM"]
    # feature_df["BRCA1&BRCA2&RAD51C&FANCF"] = feature_df["BRCA1&BRCA2&RAD51C"] + feature_df["FANCF"]
    feature_df = feature_df.sort_index()
    feature_df.to_csv(snakemake.output[0], sep="\t")
