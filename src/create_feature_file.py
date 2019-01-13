import pandas as pd
import yaml

if __name__ == '__main__':
    mc_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    mc_df.index = [i[0:12] for i in mc_df.index]
    # drop this sample because we suspect it is a melanoma metastasis instead of LUSC primary tumor
    print(mc_df.shape)
    mc_df = mc_df.drop(index=["TCGA-18-3409", "TCGA-90-A4ED"], errors="ignore")
    print(mc_df.shape)
    feature_dfs = []
    for i in range(1, len(snakemake.input)):
        feature_dfs.append(pd.read_csv(snakemake.input[i], sep="\t", index_col=0))
    # add cancer type as a covariate
    feature_dfs.append(pd.DataFrame(data=snakemake.wildcards["cancer_type"], index=mc_df.index, columns="CancerType"))
    # add number of mutations as a covariate
    feature_dfs.append(mc_df.sum(axis=1).to_frame(name="nMuts"))
    feature_df = pd.concat(feature_dfs, axis=1, join="inner")
    samples = list(set(mc_df.index).intersection(feature_df.index))
    mc_df = mc_df.loc[samples]
    mc_df = mc_df.sort_index()
    mc_df.to_csv(snakemake.output[0], sep="\t")
    feature_df = feature_df.loc[samples]
    feature_df = feature_df.sort_index()
    feature_df.to_csv(snakemake.output[1], sep="\t")
