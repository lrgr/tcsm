import pandas as pd
import yaml

if __name__ == '__main__':
    mc_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    mc_df.index = [i[0:12] for i in mc_df.index]
    # load the cancer hierarchy features
    with open(snakemake.input[1], 'r') as stream:
        cancer_hierarchy = yaml.load(stream)
    # print(cancer_hierarchy[snakemake.wildcards["cancer_type"]])
    feature_dfs = []
    for i in range(2, len(snakemake.input)):
        feature_dfs.append(pd.read_csv(snakemake.input[i], sep="\t", index_col=0))
    feature_dfs.append(pd.DataFrame(data=1, index=mc_df.index, columns=cancer_hierarchy[snakemake.wildcards["cancer_type"]]))
    feature_df = pd.concat(feature_dfs, axis=1, join="inner")
    samples = list(set(mc_df.index).intersection(feature_df.index))
    mc_df = mc_df.loc[samples]
    mc_df = mc_df.sort_index()
    mc_df.to_csv(snakemake.output[0], sep="\t")
    feature_df = feature_df.loc[samples]
    feature_df = feature_df.sort_index()
    feature_df.to_csv(snakemake.output[1], sep="\t")
