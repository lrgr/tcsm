import pandas as pd


if __name__ == '__main__':
    lst_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    lst_df.to_csv(snakemake.output[0], sep="\t")
