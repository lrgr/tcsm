import pandas as pd


if __name__ == '__main__':
    lst_df = pd.read_csv(snakemake.input[0], sep="\t", header=None, names=["Cancer Type", "Sample", "LST"])
    lst_df = lst_df.drop(columns="Cancer Type")
    lst_df = lst_df.set_index("Sample")
    lst_df.to_csv(snakemake.output[0], sep="\t")
