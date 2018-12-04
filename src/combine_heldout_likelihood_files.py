import pandas as pd
import re

if __name__ == '__main__':
    #print(snakemake.wildcards["n_samples"])
    #print(len(snakemake.input))
    output = []
    for i in range(0, len(snakemake.input)):
        df = pd.read_csv(snakemake.input[i], sep="\t")
        print(df)
        output.append(df.iloc[0])
    df = pd.concat(output, axis=1)
    # print(df)
    df = df.transpose()
    df = df.astype({"likelihood": "float32", "K": "int32"})
    print(df)
    df = df.groupby(by="K").agg("mean")
    print(df)
    #df = df.set_index(keys="K", drop=True)
    #print(df)
    df.to_csv(snakemake.output[0], sep="\t")
    # print(df)
    # selected_K = df.loc[df["log-likelihood"].idxmax()]
    # selected_K.to_frame().to_csv(snakemake.output[1], sep="\t", index=False)
