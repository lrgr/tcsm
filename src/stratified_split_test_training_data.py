import pandas as pd
from sklearn.model_selection import StratifiedKFold
import itertools

if __name__ == '__main__':
    # first, get the number of folds
    # then, use this information to calculate the size of the test set
    mc_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    feature_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    random_state = int(snakemake.params["seed"])
    n_folds = int(snakemake.params["n_folds"])
    fold = int(snakemake.wildcards["fold"])
    sss = StratifiedKFold(n_splits=n_folds, random_state=random_state)
    all_splits = sss.split(feature_df.index, feature_df[snakemake.params["feature"]])
    train_index, test_index = next(itertools.islice(all_splits, fold, fold+1))
    mc_df.iloc[test_index].to_csv(snakemake.output[0], sep="\t")
    mc_df.iloc[train_index].to_csv(snakemake.output[1], sep="\t")
    feature_df.iloc[test_index].to_csv(snakemake.output[2], sep="\t")
    feature_df.iloc[train_index].to_csv(snakemake.output[3], sep="\t")
