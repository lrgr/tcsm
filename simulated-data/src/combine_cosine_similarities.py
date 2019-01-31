import pandas as pd
import re

if __name__ == '__main__':
    output = []
    for f in snakemake.input:
        # extract the expected number of mutations from the filename
        project = re.split("_|\.", f)[3]
        print(project)
        n_samples = re.split("-|\.", project)[0]
        iter_num = re.split("-|\.", project)[2]
        df = pd.read_csv(f, sep="\t", index_col=0)
        # get the max value
        my_df = df.max(axis=0).to_frame()
        my_df.columns = ["cosine similarity"]
        my_df["n_samples"] = n_samples
        my_df["iter_num"] = iter_num
        my_df["Signature"] = my_df.index
        output.append(my_df)

    df = pd.concat(output)
    df = df.groupby(by=["Signature", "n_samples"]).agg({"cosine similarity": 'mean'})
    df = df.unstack()
    df.columns = df.columns.droplevel(level=0)
    df.to_csv(snakemake.output[0], sep="\t")
