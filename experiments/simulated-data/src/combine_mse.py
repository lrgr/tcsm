import pandas as pd
import numpy as np
import re
from sklearn.metrics import mean_squared_error


if __name__ == '__main__':
    output = []
    for i in range(0, len(snakemake.input)):
        # extract the expected number of mutations from the filename
        n_samples = re.split("_|\.", snakemake.input[i])[2]
        iter_num = re.split("_|\.", snakemake.input[i])[4]
        # print(n_mutations)
        # print(snakemake.input[i])
        df = pd.read_csv(snakemake.input[i], sep="\t", index_col=0)
        # print(df)
        # get the max value
        df = df.apply(lambda x: mean_squared_error(x, pd.Series(0, index=x.index)), axis=0).to_frame()
        df.columns = ["mse"]
        df["n_samples"] = n_samples
        df["iter_num"] = iter_num
        df["Signature"] = df.index
        output.append(df)
    df = pd.concat(output)
    df = df.groupby(by=["Signature", "n_samples"]).agg({"mse": 'mean'})
    df = df.unstack()
    df.columns = df.columns.droplevel(level=0)
    df.to_csv(snakemake.output[0], sep="\t")
