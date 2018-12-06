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
MUTATION_COUNT_MATRIX_FILE=join(OUTPUT_DIR, "{project}_mutation_count.tsv")
FEATURE_FILE=join(OUTPUT_DIR, "{project}_features.tsv")
TRAIN_MC_FILE=join(OUTPUT_DIR, "mutation-count-train_{project}.tsv")
TRAIN_FEATURE_FILE=join(OUTPUT_DIR, "features-train_{project}.tsv")
TEST_MC_FILE=join(OUTPUT_DIR, "mutation-count-test_{project}.tsv")
TEST_FEATURE_FILE=join(OUTPUT_DIR, "features-test_{project}.tsv")

# Output files

STM_NORMALIZED_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-normalized-exposures_{covariates}_{K}_{project}.tsv')
STM_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-exposures_{covariates}_{K}_{project}.tsv')
STM_SIGNATURES_FILE=join(OUTPUT_DIR, 'stm-signatures_{covariates}_{K}_{project}.tsv')
STM_HELDOUT_LIKELIHOOD_FILE=join(OUTPUT_DIR, 'stm-heldout-likelihood_{covariates}_{K}_{project}.tsv')
# STM_MODEL_SELECTION_FILE=join(OUTPUT_DIR, 'stm-model-selection_{project}.pdf')
# STM_PERMUTATION_TEST_FILE=join(OUTPUT_DIR, 'stm-permutation-test-results_{K}_{project}.pdf')
STM_HELDOUT_LIKELIHOOD_RATIO_FILE=join(OUTPUT_DIR, 'stm-heldout-ratio_{covariates}_{K}_{project}.tsv')
STM_TEST_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-test-exposures_{covariates}_{K}_{project}.tsv')
STM_TRAIN_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-train-exposures_{covariates}_{K}_{project}.tsv')

COMBINED_HELDOUT_LIKELIHOOD_FILE=join(OUTPUT_DIR, 'stm-combined-heldout_{covariates}_{project}.tsv')
LIKELIHOOD_PLOT=join(OUTPUT_DIR,'likelihood-plot_{covariates1}_{covariates2}_{project}.pdf')
# Parameters
seed="123456"

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
        STM_TRAIN_EXPOSURES_FILE,
        STM_TEST_EXPOSURES_FILE
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
        MUTATION_COUNT_MATRIX_FILE,
        FEATURE_FILE
    output:
        STM_NORMALIZED_EXPOSURES_FILE,
        STM_SIGNATURES_FILE
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
