from sklearn.svm import LinearSVC
from scipy import interp
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from statistics import mean
from sklearn.metrics import precision_recall_curve, average_precision_score
import matplotlib.pyplot as plt


if __name__ == '__main__':
    exposure_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    feature_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    feature_df = feature_df["BRCA12"]
    # print(feature_df)
    skf = StratifiedKFold(n_splits=5, random_state=123456)
    roc_auc_scores = []
    accuracy_scores = []
    average_precision_scores = []
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    average_precisions = []
    i = 0
    for train_index, test_index in skf.split(exposure_df, feature_df):
        X_train, X_test = exposure_df.iloc[train_index], exposure_df.iloc[test_index]
        y_train, y_test = feature_df.iloc[train_index], feature_df.iloc[test_index]

        clf = LinearSVC(random_state=0, tol=1e-5, class_weight="balanced")
        clf.fit(X_train, y_train)
        # print(clf.coef_)
        # print(clf.intercept_)
        y_scores = clf.predict(X_test)
        precision, recall, thresholds = precision_recall_curve(y_test, y_scores)
        # calculate precision-recall AUC
        auc = auc(recall, precision)
        ap = average_precision_score(y_test, y_scores)
        average_precisions.append(ap)
        plt.plot(recall, precision, lw=1, alpha=0.3, label='PR fold %d (AUC = %0.2f)' % (i, roc_auc)))
        # Compute ROC curve and area the curve
    #     fpr, tpr, thresholds = roc_curve(y_test, y_scores)
    #     tprs.append(interp(mean_fpr, fpr, tpr))
    #     tprs[-1][0] = 0.0
    #     roc_auc = auc(fpr, tpr)
    #     aucs.append(roc_auc)
    #     plt.plot(fpr, tpr, lw=1, alpha=0.3,
    #          label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
    #     i += 1
    # plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
    #      label='Chance', alpha=.8)

    print(mean(average_precisions))
    # mean_tpr = np.mean(tprs, axis=0)
    # mean_tpr[-1] = 1.0
    # mean_auc = auc(mean_fpr, mean_tpr)
    # std_auc = np.std(aucs)
    # plt.plot(mean_fpr, mean_tpr, color='b',
    #          label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
    #          lw=2, alpha=.8)
    #
    # std_tpr = np.std(tprs, axis=0)
    # tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    # tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    # plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
    #                  label=r'$\pm$ 1 std. dev.')
    #
    # plt.xlim([-0.05, 1.05])
    # plt.ylim([-0.05, 1.05])
    # plt.xlabel('False Positive Rate')
    # plt.ylabel('True Positive Rate')
    # plt.title('Receiver operating characteristic example')
    # plt.legend(loc="lower right")
    # # plt.show()
    # plt.savefig(snakemake.output[0])
        # roc_auc_scores.append(roc_auc_score(y_test, y_scores)) # roc auc
        # accuracy_scores.append(clf.score(X_test, y_test)) # accuracy
        # average_precision_scores.append(average_precision_score(y_test, y_scores))
        # print(clf.predict([[0, 0, 0, 0]]))
    # print(average_precision_scores)
    # print(mean(roc_auc_scores))
    # print(mean(accuracy_scores))
    # print(mean(average_precision_scores))

    # pd.DataFrame(data={"roc auc": roc_auc_scores, "accuracy": accuracy_scores, "precision": average_precision_scores}).to_csv(snakemake.output[0])
