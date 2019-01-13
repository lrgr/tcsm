import pandas as pd, numpy as np
import math
from scipy.stats import ranksums

def simulate_logistic_normal(mean, cov, n_draws):
    y = np.random.multivariate_normal(mean, cov, size=n_draws)
    # now figure out how to convert from a multivariate_normal distribution to a logistic_normal distribution
    y = pd.DataFrame(data=y, columns=["Topic{}".format(i) for i in range(1, np.size(mean)+1)])
    y["Topic{}".format(np.size(mean)+1)] = 0
    # first take exp(e) for every element e in the data frame
    y = y.applymap(math.exp)
    # normalize each element by the sum of its row
    return y.divide(y.sum(axis=1), axis='index')


if __name__ == '__main__':
    # read input files
    sigma_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    gamma_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    # set parameters
    n_draws = 100
    np.random.seed(int(snakemake.params["seed"]))
    # assume that there is only a single binary covariate that we are interested in
    default = 1 # the first row is the intercept (1-indexed because R)
    covariate = 2
    # mean is just the default
    mean = gamma_df.loc[default].values
    cov = sigma_df.values
    topics = ["Topic{}".format(i) for i in range(1, np.size(mean)+2)]
    default_exposures = simulate_logistic_normal(mean, cov, n_draws)
    # mean is the sum of the default plus the covariate
    mean = gamma_df.loc[default].values +  gamma_df.loc[covariate].values
    cov = sigma_df.values
    covariate_exposures = simulate_logistic_normal(mean, cov, n_draws)
    for topic in topics:
        print(topic)
        print(ranksums(default_exposures[topic], covariate_exposures[topic]))
