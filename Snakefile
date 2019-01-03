# Modules
from os.path import join, dirname


wildcard_constraints:
    K="\d+",
    iter_num="\d+"

DATA_DIR = 'raw'
OUTPUT_DIR = 'processed'
SRC_DIR = 'src'

FUNCTION_FILE = join(dirname(srcdir("Snakefile")), SRC_DIR, "heldout_functions.R")

# Input Files
MC_FILE=join(OUTPUT_DIR, "mutation-count-{partition}_{project}.tsv")
ALL_MC_FILE=join(OUTPUT_DIR, "mutation-count-all_{project}.tsv")
FEATURE_FILE=join(OUTPUT_DIR, "features-all_{project}.tsv")
TRAIN_MC_FILE=join(OUTPUT_DIR, "mutation-count-train_{project}.tsv")
TRAIN_FEATURE_FILE=join(OUTPUT_DIR, "features-train_{project}.tsv")
TEST_MC_FILE=join(OUTPUT_DIR, "mutation-count-test_{project}.tsv")
TEST_FEATURE_FILE=join(OUTPUT_DIR, "features-test_{project}.tsv")
EXOME_TRINUCLEOTIDE_COUNTS = join(dirname(srcdir("Snakefile")), "exome_trinucleotide_counts.tsv")
GENOME_TRINUCLEOTIDE_COUNTS = join(dirname(srcdir("Snakefile")), "genome_trinucleotide_counts.tsv")

# Output files
# Exposures and Signatures
STM_NORMALIZED_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-{partition}-normalized-exposures_{covariates}_{K}_{project}.tsv')
STM_COUNT_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-{partition}-count-exposures_{covariates}_{K}_{project}.tsv')
STM_ALL_NORMALIZED_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-all-normalized-exposures_{covariates}_{K}_{project}.tsv')
STM_TEST_NORMALIZED_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-test-normalized-exposures_{covariates}_{K}_{project}.tsv')
STM_TRAIN_NORMALIZED_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-train-normalized-exposures_{covariates}_{K}_{project}.tsv')
STM_TRAIN_SIGNATURES_FILE=join(OUTPUT_DIR, 'stm-train-signatures_{covariates}_{K}_{project}.tsv')
STM_ALL_COUNT_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-all-count-exposures_{covariates}_{K}_{project}.tsv')
STM_TEST_COUNT_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-test-count-exposures_{covariates}_{K}_{project}.tsv')
STM_TRAIN_COUNT_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-train-count-exposures_{covariates}_{K}_{project}.tsv')

STM_EXOME_SIGNATURES_FILE=join(OUTPUT_DIR, 'stm-exome-signatures_{covariates}_{K}_{project}.tsv')
STM_GENOME_SIGNATURES_FILE=join(OUTPUT_DIR, 'stm-genome-signatures_{covariates}_{K}_{project}.tsv')
# Likelihood files
STM_HELDOUT_LIKELIHOOD_FILE=join(OUTPUT_DIR, 'stm-heldout-likelihood_{covariates}_{K}_{project}.tsv')
STM_EFFECT_TABLE=join(OUTPUT_DIR, 'stm-effect-table_{covariates}_{K}_{project}.tsv')
STM_EFFECT_PLOT=join(OUTPUT_DIR, 'stm-effect-table_{covariates}_{K}_{project}.png')
STM_HELDOUT_LIKELIHOOD_RATIO_FILE=join(OUTPUT_DIR, 'stm-heldout-ratio_{covariates}_{K}_{project}.tsv')
COMBINED_HELDOUT_LIKELIHOOD_FILE=join(OUTPUT_DIR, 'stm-combined-heldout_{covariates}_{project}.tsv')
LIKELIHOOD_PLOT=join(OUTPUT_DIR,'likelihood-plot_{covariates1}_{covariates2}_{project}.pdf')
# Parameters
seed="123456"

rule convert_normalized_exposures_to_counts:
    input:
        STM_NORMALIZED_EXPOSURES_FILE,
        MC_FILE
    output:
        STM_COUNT_EXPOSURES_FILE
    script:
        "src/convert_normalized_to_count_exposures.py"

rule convert_exome_sigs_to_genome_sigs:
    input:
        EXOME_TRINUCLEOTIDE_COUNTS,
        GENOME_TRINUCLEOTIDE_COUNTS,
        STM_EXOME_SIGNATURES_FILE
    output:
        STM_GENOME_SIGNATURES_FILE
    script:
        "src/normalize_exome_sigs_to_genome_sigs.py"

rule plot_likelihoods:
    input:
        COMBINED_HELDOUT_LIKELIHOOD_FILE.format(covariates="{covariates1}", project="{project}"),
        COMBINED_HELDOUT_LIKELIHOOD_FILE.format(covariates="{covariates2}", project="{project}")
    output:
        LIKELIHOOD_PLOT
    script:
        "src/plot_likelihoods.py"

rule stm_heldout_exposures:
    params:
        seed,
        FUNCTION_FILE
    input:
        TRAIN_MC_FILE,
        TEST_MC_FILE,
        TRAIN_FEATURE_FILE,
        TEST_FEATURE_FILE
    output:
        STM_TRAIN_NORMALIZED_EXPOSURES_FILE,
        STM_TEST_NORMALIZED_EXPOSURES_FILE,
        STM_TRAIN_SIGNATURES_FILE
    script:
        "src/stm_heldout_exposures.R"

rule stm_heldout_likelihood_ratio:
    params:
        seed,
        FUNCTION_FILE
    input:
        TRAIN_MC_FILE,
        TEST_MC_FILE,
        TRAIN_FEATURE_FILE,
        TEST_FEATURE_FILE
    output:
        STM_HELDOUT_LIKELIHOOD_RATIO_FILE
    script:
        "src/stm_heldout_likelihood_ratio.R"


rule run_stm:
    params:
        seed
    input:
        ALL_MC_FILE,
        FEATURE_FILE
    output:
        STM_ALL_NORMALIZED_EXPOSURES_FILE,
        STM_EXOME_SIGNATURES_FILE,
        STM_EFFECT_TABLE
    script:
        "src/run_stm.R"


rule stm_heldout_likelihood:
    params:
        seed,
        FUNCTION_FILE
    input:
        TRAIN_MC_FILE,
        TEST_MC_FILE,
        TRAIN_FEATURE_FILE,
        TEST_FEATURE_FILE
    output:
        STM_HELDOUT_LIKELIHOOD_FILE
    script:
        "src/stm_heldout_likelihood.R"

rule stm_plot_effect:
    input:
        STM_EFFECT_TABLE
    output:
        STM_EFFECT_PLOT
    script:
        "src/plot_effects.py"
# rule stm_permutation_test:
#     params:
#         seed,
#         COVARIATE_FORMULA
#     input:
#         MUTATION_COUNT_MATRIX_FILE,
#         FEATURE_FILE
#     output:
#         STM_PERMUTATION_TEST_FILE
#     script:
#         "src/stm_permutation_test.R"
