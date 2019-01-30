import pandas as pd
from sklearn.model_selection import StratifiedShuffleSplit

if __name__ == '__main__':
    mc_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    feature_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    random_state = int(snakemake.params["seed"])
    sss = StratifiedShuffleSplit(n_splits=1, test_size=.25, random_state=random_state)
    assert((mc_df.index == feature_df.index).all())
    for train_index, test_index in sss.split(feature_df.index, feature_df[snakemake.params["feature"]]):
        mc_df.iloc[test_index].to_csv(snakemake.output[0], sep="\t")
        mc_df.iloc[train_index].to_csv(snakemake.output[1], sep="\t")
        feature_df.iloc[test_index].to_csv(snakemake.output[2], sep="\t")
        feature_df.iloc[train_index].to_csv(snakemake.output[3], sep="\t")
