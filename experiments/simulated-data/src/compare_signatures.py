import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity


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
    df.to_csv(snakemake.output[0], sep="\t")
    # reorder exposures based on signature cosine similarity
    # in the signatures file, signatures are named "V1", "V2", "V3", "V4"
    # in the exposures file, signatures are named "Topic1", "Topic2", "Topic3", "Topic4"
    # every signature is mapped to the COSMIC signature with the highest cosine similarity
    # idx = df.idxmax(axis=1) # takes the max value for each row
    # exposures = pd.read_csv(snakemake.input[2], sep="\t", index_col=0)
    # exposures.columns = [idx.loc[col] for col in exposures.columns]
    # exposures.to_csv(snakemake.output[1], sep="\t")
