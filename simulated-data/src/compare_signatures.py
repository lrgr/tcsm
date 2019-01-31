import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity


if __name__ == '__main__':
    exp_sigs = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    gt_sigs = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    # cosine_similarity expects matrices to be samples-by-features so signatures-by-categories
    exp_sigs = exp_sigs.reindex(columns=gt_sigs.columns)
    exp_sigs = exp_sigs.fillna(0)
    M = cosine_similarity(exp_sigs, gt_sigs)
    df = pd.DataFrame(data=M, index=exp_sigs.index, columns=gt_sigs.index)
    df.to_csv(snakemake.output[0], sep="\t")
