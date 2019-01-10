import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from statistics import mean
import matplotlib.pyplot as plt
from linearSVCPR import pr_linear_SVC


if __name__ == '__main__':
    n_folds = int(len(snakemake.input)/4)
    pr_aucs = []
    for fold in range(0, n_folds):
        X_train = pd.read_csv(snakemake.input[fold], sep="\t", index_col=0)
        X_test = pd.read_csv(snakemake.input[1*n_folds+fold], sep="\t", index_col=0)
        y_train = pd.read_csv(snakemake.input[2*n_folds+fold], sep="\t", index_col=0)["BRCA1BRCA2RAD51C"]
        y_test = pd.read_csv(snakemake.input[3*n_folds+fold], sep="\t", index_col=0)["BRCA1BRCA2RAD51C"]
        pr_aucs.append(pr_linear_SVC(X_train, X_test, y_train, y_test))
        # plt.plot(recall, precision, lw=1, alpha=0.3, label='PR fold %d (AUC = %0.2f)' % (fold, pr_auc))

    print(mean(pr_aucs))
    df = pd.DataFrame(data=pr_aucs, index=range(0, n_folds), columns=["pr auc"])
    df.to_csv(snakemake.output[0], sep="\t")
    #plt.savefig(snakemake.output[0])
