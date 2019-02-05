import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity


if __name__ == '__main__':
    # first, read in the estimated exposures to estimated signatures
    est_exposures = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    # next, read in the estimated signatures
    est_sigs = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    # next, read in the ground truth signatures
    gt_sigs = pd.read_csv(snakemake.input[2], sep="\t", index_col=0)
    # read in the ground truth exposures
    gt_exposures = pd.read_csv(snakemake.input[3], sep="\t", index_col=0)
    # now, we need to align the columns between the estimated and ground truth exposures
    est_sigs = est_sigs.reindex(columns=gt_sigs.columns)
    est_sigs = est_sigs.fillna(0)
    M = cosine_similarity(est_sigs, gt_sigs)
    df = pd.DataFrame(data=M, index=est_sigs.index, columns=gt_sigs.index)
    #print(df.max(axis=0)) # max for each column
    #print(df.idxmax(axis=0)) # returns the topic that most resembles each ground truth signature
    #print(df.max(axis=1)) # max for each row
    idx = df.idxmax(axis=1) # returns the ground truth signature that most resembles each topic
    # corner case when multiple extracted signatures have the highest cosine similarity to one ground truth signature
    if (len(idx) != len(set(idx))):
        # first, we need to identify the extracted signatures that map to the same ground truth signature
        duplicated_topics = idx.duplicated(keep=False)
        # now get the cosine similarity dataframe for just the duplicated topics
        duplicated_topics_df = df.loc[duplicated_topics]
        # get the topic with the lower max cosine similarity
        min_topic = duplicated_topics_df.max(axis=1).idxmin()
        # get the ground truth signature that wasn't assigned
        missing_gt_signature = set(gt_sigs.index).difference(set(idx))
        missing_gt_signature = list(missing_gt_signature)[0]
        #print(min_topic)
        #print(missing_gt_signature)
        idx[min_topic] = missing_gt_signature

        # between the two extracted signatures, we want to assign the ground signature to the one that matches best

        # we want to assign the signature with the

    assert(len(idx) == len(set(idx)))
    # idx is a series where the signatures are ordered correctly
    est_exposures.columns = idx.values #[idx.loc[col] for col in est_exposures.columns]
    diff_exposures = gt_exposures.subtract(est_exposures)
    n_mutations = gt_exposures.sum(axis='columns')
    # print(n_mutations)
    norm_diff_exposures = diff_exposures.div(n_mutations, axis='index')
    # print(diff_exposures)
    # print(norm_diff_exposures)
    norm_diff_exposures.to_csv(snakemake.output[0], sep="\t")
