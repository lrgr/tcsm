import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib as plt
import seaborn as sns


if __name__ == '__main__':
    exp_sigs = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    # exp_sigs = exp_sigs.transpose()
    #exp_sigs = exp_sigs.sort_index(axis=1)
    # exp_sigs.index = [i.replace("V", "Topic") for i in exp_sigs.index]
    gt_sigs = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    #gt_sigs = gt_sigs.sort_index(axis=1)

    # cosine_similarity expects matrices to be samples-by-features so signatures-by-categories
    exp_sigs = exp_sigs.reindex(columns=gt_sigs.columns)
    #df = pd.concat([exp_sigs, gt_sigs])
    exp_sigs = exp_sigs.fillna(0)
    M = cosine_similarity(exp_sigs, gt_sigs)
    df = pd.DataFrame(data=M, index=exp_sigs.index, columns=gt_sigs.index)
    df.to_csv(snakemake.output[1], sep="\t")
    ax = sns.heatmap(df, annot=True)
    fig = ax.get_figure()
    fig.savefig(snakemake.output[0])
    # print(df)
    # df.to_csv(snakemake.output[0], sep="\t")
