import pandas as pd
from statsmodels.distributions.empirical_distribution import ECDF

if __name__ == '__main__':
    cov_sig = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    shuffled_cov_sig_list = []
    for i in range(1, len(snakemake.input)):
        df = pd.read_csv(snakemake.input[i], sep="\t", index_col=0)
        shuffled_cov_sig_list.append(df["covariate mean"] - df["default mean"])
    df = pd.concat(shuffled_cov_sig_list, axis=1)
    cov_sig = cov_sig["covariate mean"] - cov_sig["default mean"]
    output = []
    for topic in df.index:
        cdf = ECDF(df.loc[topic])
        output.append(cdf(cov_sig.loc[topic]))
    pd.Series(index=df.index, data=output).to_csv(snakemake.output[0], sep="\t")
