import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import string
import sys

sys.path.append(snakemake.params[0])
from plot_likelihoods import plot_likelihoods, combine_likelihood_files
from plot_heldout_ratio import combine_ratio_files, plot_heldout_ratio
from plot_LST import plot_LST
fig, axs = plt.subplots(1,3,figsize=(18,6))

# first, plot the heldout log-likelihood across K
df = combine_likelihood_files(snakemake.input["likelihood_files"])
plot_likelihoods(df, axs[0])
axs[0].text(-0.1, 1.1, "A", transform=axs[0].transAxes,
            size=20, weight='bold')

# second, plot the heldout log-likelihood ratio for the HR covariate
df = combine_ratio_files(snakemake.input["ratio_files"])
plot_heldout_ratio(df, axs[1])
axs[1].text(-0.1, 1.1, "B", transform=axs[1].transAxes,
            size=20, weight='bold')

# third, plot the LST count vs heldout log-likelihood ratio
axs[2].text(-0.1, 1.1, "C", transform=axs[2].transAxes,
            size=20, weight='bold')
plot_LST(df, axs[2])
plt.savefig(snakemake.output[0])
plt.savefig(snakemake.output[1])
