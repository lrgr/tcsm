import pandas as pd
import numpy as np
import re
from sklearn.metrics import mean_squared_error


if __name__ == '__main__':
    output = []
    for f in snakemake.input:
        # extract the expected number of mutations from the filename
        project = re.split("_|\.", f)[3]
        n_samples = re.split("-|\.", project)[0]
        iter_num = re.split("-|\.", project)[2]
        df = pd.read_csv(f, sep="\t", index_col=0)
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
