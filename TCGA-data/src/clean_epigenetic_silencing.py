import pandas as pd

if __name__ == '__main__':
    df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    df.columns = ["es{}".format(g) for g in df.columns]
    df.to_csv(snakemake.output[0], sep="\t")
