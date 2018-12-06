import pandas as pd

if __name__ == '__main__':
    mc_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    clin_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)

    lst_df = pd.read_csv(snakemake.input[2], sep="\t", header=None, names=["Cancer Type", "Sample", "LST"])
    lst_df = lst_df.drop(columns="Cancer Type")
    lst_df = lst_df.set_index("Sample")
    mc_df.index = [i[0:12] for i in mc_df.index]
    samples = list(set(mc_df.index).intersection(clin_df.index).intersection(lst_df.index))
    mc_df.loc[samples].to_csv(snakemake.output[0], sep="\t")
    clin_df = clin_df.loc[samples]
    lst_df = lst_df.loc[samples]
    clin_df = clin_df.merge(lst_df, left_index=True, right_index=True)
    clin_df.to_csv(snakemake.output[1], sep="\t")
