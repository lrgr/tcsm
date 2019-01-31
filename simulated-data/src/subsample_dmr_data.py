import pandas as pd
# import numpy as np


if __name__ == '__main__':
    f_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    mc_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    e_df = pd.read_csv(snakemake.input[2], sep="\t", index_col=0)
    n = int(snakemake.wildcards["n_samples"])
    f_df.iloc[0:n].to_csv(snakemake.output[0], sep="\t")
    mc_df.iloc[0:n].to_csv(snakemake.output[1], sep="\t")
    e_df.iloc[0:n].to_csv(snakemake.output[2], sep="\t")
