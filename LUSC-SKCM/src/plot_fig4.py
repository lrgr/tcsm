import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import string
import sys

sys.path.append(snakemake.params[0])
from plot_likelihoods import plot_likelihoods, combine_likelihood_files
from plot_ratio_dnp import plot_likelihood_ratio

fig, axs = plt.subplots(1,2,figsize=(12,6))

# first, plot the heldout log-likelihood across K
df = combine_likelihood_files(snakemake.input[1:])

plot_likelihoods(df, axs[0])
axs[0].text(-0.1, 1.1, "A", transform=axs[0].transAxes,
            size=20, weight='bold')
# second, plot the heldout log-likelihood for the two cancer types
df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
plot_likelihood_ratio(df, axs[1])
axs[1].text(-0.1, 1.1, "B", transform=axs[1].transAxes,
            size=20, weight='bold')
plt.savefig(snakemake.output[0])
plt.savefig(snakemake.output[1])
