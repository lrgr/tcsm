import pandas as pd

if __name__ == '__main__':
    bi_inactivation_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    bi_inactivation_df.columns = ["bi{}".format(g) for g in bi_inactivation_df.columns]
    bi_inactivation_df.to_csv(snakemake.output[0], sep="\t")
