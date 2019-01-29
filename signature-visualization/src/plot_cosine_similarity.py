import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib as plt
import seaborn as sns; sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
# rcParams['font.family'] = 'Courier New'


if __name__ == '__main__':
    exp_sigs = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    gt_sigs = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    # cosine_similarity expects matrices to be samples-by-features so signatures-by-categories
    exp_sigs = exp_sigs.reindex(columns=gt_sigs.columns)
    #df = pd.concat([exp_sigs, gt_sigs])
    exp_sigs = exp_sigs.fillna(0)
    M = cosine_similarity(exp_sigs, gt_sigs)
    df = pd.DataFrame(data=M, index=range(1, len(exp_sigs.index)+1), columns=range(1, len(gt_sigs.index)+1))
    df.to_csv(snakemake.output[1], sep="\t")
    ax = sns.heatmap(df, annot=True, cbar_kws={'label': 'Cosine similarity'})
    ax.set_ylabel("COSMIC")
    ax.set_xlabel("TCSM", fontname='Courier New')
    fig = ax.get_figure()

    fig.savefig(snakemake.output[0])
