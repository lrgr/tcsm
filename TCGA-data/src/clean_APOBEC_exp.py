import pandas as pd

if __name__ == '__main__':
    apobec_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    apobec_df = apobec_df.drop(columns=["gene_id.1"])
    apobec_df = apobec_df.transpose()
    apobec_df.index = [s[0:15] for s in apobec_df.index]
    mc_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    mc_df.index = [s[0:15] for s in mc_df.index]
    apobec_df = apobec_df.loc[apobec_df.index.isin(mc_df.index)]
    # print(apobec_df.index.unique().shape)
    apobec_df.index = [s[0:12] for s in apobec_df.index]
    apobec_df = apobec_df.drop(index=["TCGA-18-3409", "TCGA-90-A4ED"], errors="ignore")
    # print(apobec_df.index.drop_duplicates())
    apobec_df = apobec_df[~apobec_df.index.duplicated()]
    # print(apobec_df.index.shape)
    # print(apobec_df.index.unique().shape)
    apobec_df.to_csv(snakemake.output[0], sep="\t")
