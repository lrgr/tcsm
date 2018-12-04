# Modules
from os.path import join, dirname


wildcard_constraints:
    K="\d+",
    iter_num="\d+",
    model="stm|ctm"

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
HELDOUT_LIKELIHOOD_FILE=join(OUTPUT_DIR, '{model}-heldout-likelihood_{K}_{project}.tsv')


# CTM_NORMALIZED_EXPOSURES_FILE=join(OUTPUT_DIR, 'ctm-normalized-exposures_{K}_{project}.tsv')
# CTM_EXPOSURES_FILE=join(OUTPUT_DIR, 'ctm-exposures_{K}_{project}.tsv')
# CTM_SIGNATURES_FILE=join(OUTPUT_DIR, 'ctm-signatures_{K}_{project}.tsv')
# CTM_HELDOUT_LIKELIHOOD_FILE=join(OUTPUT_DIR, 'ctm-heldout-likelihood_{K}_{project}.tsv')
# CTM_HELDOUT_BOUND_FILE=join(OUTPUT_DIR, 'ctm-heldout-bound_{K}_{project}.tsv')
# CTM_MODEL_SELECTION_FILE=join(OUTPUT_DIR, 'ctm-model-selection_{project}.pdf')
# CTM_TEST_EXPOSURES_FILE=join(OUTPUT_DIR, 'ctm-test-exposures_{K}_{project}.tsv')
# CTM_TRAIN_EXPOSURES_FILE=join(OUTPUT_DIR, 'ctm-train-exposures_{K}_{project}.tsv')


STM_NORMALIZED_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-normalized-exposures_{covariates}_{K}_{project}.tsv')
STM_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-exposures_{covariates}_{K}_{project}.tsv')
STM_SIGNATURES_FILE=join(OUTPUT_DIR, 'stm-signatures_{covariates}_{K}_{project}.tsv')
# STM_HELDOUT_LIKELIHOOD_FILE=join(OUTPUT_DIR, 'stm-heldout-likelihood_{K}_{project}.tsv')
STM_HELDOUT_LIKELIHOOD_FILE=join(OUTPUT_DIR, 'stm-heldout-likelihood_{covariates}_{K}_{project}.tsv')
STM_HELDOUT_BOUND_FILE=join(OUTPUT_DIR, 'stm-heldout-bound_{K}_{project}.tsv')
STM_MODEL_SELECTION_FILE=join(OUTPUT_DIR, 'stm-model-selection_{project}.pdf')
STM_PERMUTATION_TEST_FILE=join(OUTPUT_DIR, 'stm-permutation-test-results_{K}_{project}.pdf')
STM_HELDOUT_LIKELIHOOD_RATIO_FILE=join(OUTPUT_DIR, 'stm-heldout-ratio_{covariates}_{K}_{project}.tsv')
STM_TEST_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-test-exposures_{K}_{project}.tsv')
STM_TRAIN_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-train-exposures_{K}_{project}.tsv')

COMBINED_HELDOUT_LIKELIHOOD_FILE=join(OUTPUT_DIR, 'combined-heldout_{model}_{project}.tsv')
LIKELIHOOD_PLOT=join(OUTPUT_DIR,'likelihood-plot_{model1}_{model2}_{project}.pdf')
# Parameters
seed="123456"
# COVARIATE_FORMULA="~ feature1"
# min_K=4
# max_K=5

rule plot_likelihoods:
    input:
        expand(COMBINED_HELDOUT_LIKELIHOOD_FILE, model="{model1}", project="{project}"),
        expand(COMBINED_HELDOUT_LIKELIHOOD_FILE, model="{model2}", project="{project}")
    output:
        LIKELIHOOD_PLOT
    script:
        "src/plot_likelihoods.py"

rule stm_heldout_exposures:
    params:
        seed,
        COVARIATE_FORMULA,
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

rule ctm_heldout_exposures:
    params:
        seed,
        FUNCTION_FILE
    input:
        TRAIN_MC_FILE,
        TEST_MC_FILE
    output:
        CTM_TRAIN_EXPOSURES_FILE,
        CTM_TEST_EXPOSURES_FILE
    script:
        "src/ctm_heldout_exposures.R"

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

rule run_ctm:
    params:
        seed
    input:
        MUTATION_COUNT_MATRIX_FILE
    output:
        CTM_NORMALIZED_EXPOSURES_FILE,
        CTM_SIGNATURES_FILE
    script:
        "src/run_ctm.R"


rule ctm_heldout_likelihood:
    params:
        seed,
        FUNCTION_FILE
    input:
        TRAIN_MC_FILE,
        TEST_MC_FILE
    output:
        CTM_HELDOUT_LIKELIHOOD_FILE
    script:
        "src/ctm_heldout_likelihood.R"

rule run_stm:
    params:
        seed,
        COVARIATE_FORMULA
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

rule stm_permutation_test:
    params:
        seed,
        COVARIATE_FORMULA
    input:
        MUTATION_COUNT_MATRIX_FILE,
        FEATURE_FILE
    output:
        STM_PERMUTATION_TEST_FILE
    script:
        "src/stm_permutation_test.R"

rule run_ctm_model_selection:
    params:
        seed
    input:
        MUTATION_COUNT_MATRIX_FILE
    output:
        CTM_MODEL_SELECTION_FILE
    script:
        "src/select_model_ctm.R"

rule run_stm_model_selection:
    params:
        seed,
        COVARIATE_FORMULA
    input:
        MUTATION_COUNT_MATRIX_FILE,
        FEATURE_FILE
    output:
        STM_MODEL_SELECTION_FILE
    script:
        "src/select_model_stm.R"

rule stm_convert_normalized_to_count_exposures:
    input:
        STM_NORMALIZED_EXPOSURES_FILE,
        MUTATION_COUNT_MATRIX_FILE
    output:
        STM_EXPOSURES_FILE
    script:
        "src/convert_normalized_to_count_exposures.py"

rule ctm_convert_normalized_to_count_exposures:
    input:
        CTM_NORMALIZED_EXPOSURES_FILE,
        MUTATION_COUNT_MATRIX_FILE
    output:
        CTM_EXPOSURES_FILE
    script:
        "src/convert_normalized_to_count_exposures.py"

rule ctm_heldout_bound:
    params:
        seed
    input:
        TRAIN_MC_FILE,
        TEST_MC_FILE
    output:
        CTM_HELDOUT_BOUND_FILE
    script:
        "src/ctm_heldout_bound.R"

rule stm_heldout_bound:
    params:
        seed,
        COVARIATE_FORMULA,
        FUNCTION_FILE
    input:
        TRAIN_MC_FILE,
        TRAIN_FEATURE_FILE,
        TEST_MC_FILE,
        TEST_FEATURE_FILE
    output:
        STM_HELDOUT_BOUND_FILE
    script:
        "src/stm_heldout_bound.R"
