import pandas as pd
import numpy as np

if __name__ == '__main__':
    df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    mc_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    # sum the rows of the mutation count matrix to get the number of mutations per sample
    n_mutations = mc_df.sum(axis=1)
    print(n_mutations)
    # df is samples-by-signatures
    # n_mutations is vector of length samples
    df = df.transpose().multiply(n_mutations).transpose()
    #df.index = [i-1 for i in df.index]
    df.to_csv(snakemake.output[0], sep="\t")
