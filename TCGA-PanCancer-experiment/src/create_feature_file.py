import pandas as pd

if __name__ == '__main__':
    mc_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    feature_df = pd.DataFrame(data=1, index=mc_df.index, columns=[snakemake.wildcards["cancer_type"]])
    feature_df.to_csv(snakemake.output[0], sep="\t")
