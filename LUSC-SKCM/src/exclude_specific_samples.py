import pandas as pd

if __name__ == '__main__':
    mc_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    feature_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    df = pd.read_excel(snakemake.input[2])
    mc_df = mc_df.loc[~mc_df.index.isin(df["Excluded"])]
    feature_df = feature_df.loc[~feature_df.index.isin(df["Excluded"])]
    if snakemake.wildcards["cancer_type"] == "SKCM":
        feature_df["Smoker"] = 0
    mc_df.to_csv(snakemake.output[0], sep="\t")
    feature_df.to_csv(snakemake.output[1], sep="\t")
