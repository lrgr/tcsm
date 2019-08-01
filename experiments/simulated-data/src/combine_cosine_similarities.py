import pandas as pd
import re

if __name__ == '__main__':
    output = []
    for i in range(0, len(snakemake.input)):
        # extract the expected number of mutations from the filename
        n_samples = re.split("_|\.", snakemake.input[i])[2]
        iter_num = re.split("_|\.", snakemake.input[i])[4]
        df = pd.read_csv(snakemake.input[i], sep="\t", index_col=0)
        # get the max value
        my_df = df.max(axis=0).to_frame()
        my_df.columns = ["cosine similarity"]
        my_df["n_samples"] = n_samples 
        my_df["iter_num"] = iter_num
        my_df["Signature"] = my_df.index
        # print(my_df)
        output.append(my_df)

    df = pd.concat(output)
    df = df.groupby(by=["Signature", "n_samples"]).agg({"cosine similarity": 'mean'})
    df = df.unstack()
    df.columns = df.columns.droplevel(level=0)
    print(df)
    print(df.columns)
    print(df.index)
    df.to_csv(snakemake.output[0], sep="\t")
