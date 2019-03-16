import pandas as pd

if __name__ == '__main__':
    seed = int(snakemake.params["seed"])+int(snakemake.wildcards["iter_num"])
    mc_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    feature_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    shuffled_feature_df = feature_df.sample(frac=1, random_state=seed)
    mc_df.to_csv(snakemake.output[0], sep="\t")
    shuffled_feature_df.to_csv(snakemake.output[1], sep="\t")
