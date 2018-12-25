#!/usr/bin/env python

# Load required modules
import matplotlib
matplotlib.use('Agg')
import sys, os, argparse, pandas as pd, numpy as np
import seaborn as sns, matplotlib.pyplot as plt
from sklearn.externals import joblib
sns.set_style('whitegrid')

# Load mutation signatures visualizations
this_dir = os.path.dirname(__file__)
viz_dir = os.path.join(this_dir, '..', '..', '../mutation-signatures-viz/src')
print(viz_dir)
sys.path.append( viz_dir )
from mutation_signatures_visualization import plot_signatures


# Parse command-line arguments
# parser = argparse.ArgumentParser()
# parser.add_argument('-cf', '--counts_file', type=str, required=True)
# parser.add_argument('-sf', '--signature_file', type=str, required=True)
# parser.add_argument('-ef', '--exposure_file', type=str, required=True)
# parser.add_argument('-of', '--output_file', type=str, required=True)
# args = parser.parse_args(sys.argv[1:])

# Load the signatures and the counts
sigs = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
# Load the exposures
exposures = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)


# Load the counts
sbs96_df = pd.read_csv(snakemake.input[2], sep='\t', index_col=0)

# required input are the counts_df, the sigs_df and the exposures_df
# plot_signatures(counts_df, signature_df, exposure_df, output_file)
plot_signatures(sbs96_df, sigs, exposures, snakemake.output[0])
# phi are the signatures
# theta is the proportional contribution of each signature to each patient
# Compute contribution per signature
