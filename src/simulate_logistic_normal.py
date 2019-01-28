import pandas as pd, numpy as np
import math
from scipy.stats import ranksums
from statsmodels.stats.multitest import fdrcorrection

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
    feature_df = pd.read_csv(snakemake.input[2], sep="\t", index_col=0)
    # set parameters
    n_draws = 5000
    # set the seed for reproducibility
    np.random.seed(int(snakemake.params["seed"]))
    # assume that there is only a single binary covariate that we are interested in
    covariate_of_interest = str(snakemake.wildcards["covariate_of_interest"])
    # identify if any other covariates exist (because if so, we need to sample them)
    covariate_list = snakemake.wildcards["covariates"].split("+")
    other_covariates = list(set(covariate_list).difference([covariate_of_interest]))
    if len(other_covariates) == 0:
        cov = sigma_df.values
        mean = gamma_df.loc["default"].values
        default_exposures = simulate_logistic_normal(mean, cov, n_draws)
        mean = gamma_df.loc["default"].values + gamma_df.loc[covariate_of_interest].values
        covariate_exposures = simulate_logistic_normal(mean, cov, n_draws)
    else:
        print(feature_df[other_covariates])
        default_features = feature_df.loc[feature_df[covariate_of_interest] == 0][other_covariates]
        # sample the feature with replacement n_draws times
        sampled_features = default_features.sample(n=n_draws, replace=True)
        mean = gamma_df.loc["default"].values
        print(mean)
        for cov in other_covariates:
            mean = mean + np.outer(sampled_features[cov].values, gamma_df.loc[cov].values)
            print(mean)
        # now we have the mean, draw from the logistic normal distribution
        output = []
        for m in mean:
            output.append(simulate_logistic_normal(mean=m, cov=sigma_df.values, n_draws=1))
        default_exposures = pd.concat(output)
        print(default_exposures)
        # for each row, calculate the mean and draw from the logistic normal
        covariate_features = feature_df.loc[feature_df[covariate_of_interest] == 1][other_covariates]
        # sample the feature with replacement n_draws times
        sampled_features = covariate_features.sample(n=n_draws, replace=True)
        mean = gamma_df.loc["default"].values + gamma_df.loc[covariate_of_interest].values
        print(mean)
        for cov in other_covariates:
            mean = mean + np.outer(sampled_features[cov].values, gamma_df.loc[cov].values)
            print(mean)
        output = []
        for m in mean:
            output.append(simulate_logistic_normal(mean=m, cov=sigma_df.values, n_draws=1))
        covariate_exposures = pd.concat(output)
        print(covariate_exposures)
    print(np.size(mean, axis=0))
    print(np.size(mean, axis=1))
    topics = ["Topic{}".format(i) for i in range(1, np.size(mean, axis=1)+2)]
    output = []
    for topic in topics:
        ranksum_result = ranksums(covariate_exposures[topic], default_exposures[topic])
        print(ranksum_result)
        data = [ranksum_result[0], ranksum_result[1], np.mean(covariate_exposures[topic]), np.mean(default_exposures[topic])]
        s = pd.Series(data=data, index=["Wilcoxon Statistic", "Wilcoxon p-value", "covariate mean", "default mean"])
        output.append(s)
    df = pd.concat(output, axis=1).transpose()
    # calculate the FDR - how?

    df["Wilcoxon FDR"] = fdrcorrection(df["Wilcoxon p-value"])[1]
    df.to_csv(snakemake.output[0], sep="\t")
