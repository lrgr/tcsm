import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from statistics import mean
from linearSVCPR import pr_linear_SVC



if __name__ == '__main__':
    exposure_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    feature_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    # exposure_df = feature_df["LST"]
    feature_df = feature_df["BRCA12"]
    skf = StratifiedKFold(n_splits=5, random_state=123456)
    pr_aucs = []
    for train_index, test_index in skf.split(exposure_df, feature_df):
        X_train, X_test = exposure_df.iloc[train_index], exposure_df.iloc[test_index]
        y_train, y_test = feature_df.iloc[train_index], feature_df.iloc[test_index]
        pr_aucs.append(pr_linear_SVC(X_train, X_test, y_train, y_test))
    print(mean(pr_aucs))
    df = pd.DataFrame(data=pr_aucs, index=range(0,5), columns=["pr auc"])
    df.to_csv(snakemake.output[0], sep="\t")
