from sklearn.svm import LinearSVC
from scipy import interp
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from statistics import mean
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, auc
import matplotlib.pyplot as plt


if __name__ == '__main__':
    n_folds = int(len(snakemake.input)/4)
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    for fold in range(0, n_folds):

        X_train = pd.read_csv(snakemake.input[fold], sep="\t", index_col=0)
        X_test = pd.read_csv(snakemake.input[1*n_folds+fold], sep="\t", index_col=0)
        y_train = pd.read_csv(snakemake.input[2*n_folds+fold], sep="\t", index_col=0)["BRCA12"]
        y_test = pd.read_csv(snakemake.input[3*n_folds+fold], sep="\t", index_col=0)["BRCA12"]
        clf = LinearSVC(random_state=0, tol=1e-5, class_weight="balanced")
        clf.fit(X_train, y_train)
        # print(clf.coef_)
        # print(clf.intercept_)
        y_scores = clf.predict(X_test)
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y_test, y_scores)
        tprs.append(interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        plt.plot(fpr, tpr, lw=1, alpha=0.3,
             label='ROC fold %d (AUC = %0.2f)' % (fold, roc_auc))

    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
         label='Chance', alpha=.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    plt.plot(mean_fpr, mean_tpr, color='b',
             label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
             lw=2, alpha=.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                     label=r'$\pm$ 1 std. dev.')

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    # plt.show()
    plt.savefig(snakemake.output[0])
    # print(feature_df)
    # clf = LinearSVC(random_state=0, tol=1e-5, class_weight="balanced")
    # clf.fit(X_train, y_train)
    # accuracy = clf.score(X_test, y_test)
    # y_scores = clf.predict(X_test)
    # roc_auc_score = roc_auc_score(y_test, y_scores) # roc auc
    # accuracy_score = clf.score(X_test, y_test) # accuracy
    # average_precision_score = average_precision_score(y_test, y_scores)
    # print(score)
    # pd.DataFrame(data={"accuracy": accuracy_score, "fold": snakemake.wildcards["fold"], "model": snakemake.wildcards["covariates"], "precision": average_precision_score, "roc auc": roc_auc_score}, index=[0]).to_csv(snakemake.output[0], sep="\t", index=False)
