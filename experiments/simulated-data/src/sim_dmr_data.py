#!/usr/bin/env python

# Load required modules
import sys, os, argparse, numpy as np, pandas as pd
from scipy.stats import multivariate_normal, poisson

if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-at', '--active_topics', type=int, required=False, default=[],
        nargs='*', help='Indices of topics to use from the topic file.') # if not specified, uses all topics
    parser.add_argument('-tf', '--topics_file', type=str, required=True)
    parser.add_argument('-n', '--n_samples', type=int, required=False, default=-1)
    parser.add_argument('-ff', '--features_file', type=str, required=True)
    parser.add_argument('-fn', '--feature_names', type=str, required=False, default=None, nargs='*')
    parser.add_argument('-l', '--lambda_file', type=str, required=True)
    parser.add_argument('-enm', '--empirical_n_mut_file', type=str, required=True)
    parser.add_argument('-om', '--output_model_file', type=str, required=True)
    parser.add_argument('-omf', '--output_mutation_count_matrix_file', type=str, required=True)
    parser.add_argument('-os', '--output_signature_file', type=str, required=True)
    parser.add_argument('-oe', '--output_exposures_file', type=str, required=True)
    parser.add_argument('-od', '--output_documents_file', type=str, required=True)
    parser.add_argument('-rs', '--random_seed', type=int, required=False, default=29810)
    parser.add_argument('-ni', '--iteration_number', type=int, required=True)
    args = parser.parse_args(sys.argv[1:])

    random_state=args.random_seed+args.iteration_number

    np.random.seed(random_state)

    # Load the given topics
    topics_df = pd.read_csv(args.topics_file, sep='\t', index_col=0)
    words = topics_df.columns
    topics = topics_df.values
    if len(args.active_topics) > 0:
        topics = topics[np.array(args.active_topics)-1]

    K, L = topics.shape
    Ks = list(range(K))

    # Load the features file
    features_df = pd.read_csv(args.features_file, sep='\t', index_col=0)
    samples = features_df.index

    if args.n_samples < len(samples):
        samples = np.random.permutation(samples)[:args.n_samples]
        features_df = features_df.loc[features_df.index.isin(set(samples))]

    N = len(samples)

    if args.feature_names is None:
        args.feature_names = features_df.columns

    features = args.feature_names
    x = features_df[features].values
    n_features = len(features)

    # Sample hyperparameters: lambda and number of mutations
    lambdas = pd.read_csv(args.lambda_file, sep="\t").values
    # sigma = 1
    # lambdas = [ multivariate_normal.rvs([0.] * n_features, np.eye(n_features) * sigma)
    #            for _ in range(K) ]
    # manually set lambas to be a
    # [1, 2, 3, 5]
    #lambdas = np.array([[.03, -4], [0, -1.5], [0, -2.5], [.035, -4]])
    df = pd.read_csv(args.empirical_n_mut_file, sep="\t", header=None, names=["n mutations"], dtype=np.int32)
    #n_mutations = poisson.rvs(args.expected_n_mut, size=N)
    n_mutations = df.sample(N, replace=True, random_state=random_state)["n mutations"].values

    # Generate data
    w = []
    z = []
    alphas = []
    thetas = []
    for i, (x_i, n_mut) in enumerate(zip(x, n_mutations)):
        # Get priors and sample mixings
        alpha_i = np.array([ np.exp(x_i.dot(lambdas_t)) for lambdas_t in lambdas ])
        theta_i = np.random.dirichlet(alpha_i)
        # print(alpha_i)
        # Sample words
        w_i = []
        z_i = []
        for j in range(n_mut):
            z_i.append(np.random.choice(Ks, p=theta_i))
            w_i.append(np.random.choice(words, p=topics[z_i[-1]]))

        # Record everything
        # print(z_i)
        w.append(w_i)
        z.append(z_i)
        thetas.append(theta_i)
        alphas.append(alpha_i)

    # Compute correlations between features and mutation signature exposure levels
    samples = [str(s) for s in samples]

    # Output in expected mutation count format
    # w is a list of lists where each list contains the mutations for a given sample
    series_list = []
    for n in range(0, len(samples)):
        series_list.append(pd.value_counts(w[n]))
    # generate a df where samples are columns and mutation categories are rows
    df = pd.concat(series_list, axis=1)
    # updates the df to have all mutation categories (even the ones that are not present)
    df = df.reindex(index=words)
    # all missing values should be 0
    df = df.fillna(0)
    # expected format is categories as columns and samples as rows
    df = df.transpose()
    df.to_csv(args.output_mutation_count_matrix_file, sep="\t")

    # Output the signatures in expected format
    signatures = ["SBS{}".format(i) for i in args.active_topics]
    signatures_df = pd.DataFrame(data=topics, columns=words, index=signatures)
    signatures_df.to_csv(args.output_signature_file, sep="\t")
    # Output to file
    with open(args.output_documents_file, 'w') as OUT:
        OUT.write('#Sample\t %s \t Words\n' % '\t'.join(features))
        for i, s in enumerate(samples):
            OUT.write(s + '\t')
            OUT.write('\t'.join(map(str, x[i])) + '\t')
            OUT.write(' '.join(w[i]) + '\n')

    sample_topic_list = []
    for n in range(0, len(samples)):
        sample_topic_list.append(pd.value_counts(z[n]))
    df = pd.concat(sample_topic_list, axis=1)
    df = df.fillna(0)
    df.index = ["SBS{}".format(i) for i in args.active_topics]
    df = df.transpose()
    df.to_csv(args.output_exposures_file, sep="\t")

    np.savez(args.output_model_file, params=vars(args),
             w=w, z=z, x=x, lambdas=lambdas, thetas=thetas, alphas=alphas,
             topics=topics, n_mutations=n_mutations)
