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
MC_FILE=join(OUTPUT_DIR, "mutation-count-all_{project}.tsv")
FEATURE_FILE=join(OUTPUT_DIR, "features-all_{project}.tsv")
TRAIN_MC_FILE=join(OUTPUT_DIR, "mutation-count-train_{project}.tsv")
TRAIN_FEATURE_FILE=join(OUTPUT_DIR, "features-train_{project}.tsv")
TEST_MC_FILE=join(OUTPUT_DIR, "mutation-count-test_{project}.tsv")
TEST_FEATURE_FILE=join(OUTPUT_DIR, "features-test_{project}.tsv")
VALIDATE_MC_FILE=join(OUTPUT_DIR, "mutation-count-validate_{project}.tsv")
VALIDATE_FEATURE_FILE=join(OUTPUT_DIR, "features-validate_{project}.tsv")

# Output files
# Exposures and Signatures
STM_NORMALIZED_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-normalized-exposures_{covariate_of_interest}_{covariates}_{K}_{project}.tsv')
STM_ALL_NORMALIZED_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-all-normalized-exposures_{covariates}_{K}_{project}.tsv')
STM_TEST_NORMALIZED_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-test-normalized-exposures_{covariate_of_interest}_{covariates}_{K}_{project}.tsv')
STM_TRAIN_NORMALIZED_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-train-normalized-exposures_{covariate_of_interest}_{covariates}_{K}_{project}.tsv')
STM_VALIDATE_NORMALIZED_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-validate-normalized-exposures_{covariate_of_interest}_{covariates}_{K}_{project}.tsv')

STM_TRAIN_SIGNATURES_FILE=join(OUTPUT_DIR, 'stm-train-signatures_{covariate_of_interest}_{covariates}_{K}_{project}.tsv')
STM_SIGNATURES_FILE=join(OUTPUT_DIR, 'stm-signatures_{covariate_of_interest}_{covariates}_{K}_{project}.tsv')

STM_ALL_COUNT_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-all-count-exposures_{covariates}_{K}_{project}.tsv')
STM_TEST_COUNT_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-test-count-exposures_{covariate_of_interest}_{covariates}_{K}_{project}.tsv')
STM_TRAIN_COUNT_EXPOSURES_FILE=join(OUTPUT_DIR, 'stm-train-count-exposures_{covariates}_{K}_{project}.tsv')

STM_EXOME_SIGNATURES_FILE=join(OUTPUT_DIR, 'stm-exome-signatures_{covariates}_{K}_{project}.tsv')

# Likelihood files
STM_HELDOUT_LIKELIHOOD_FILE=join(OUTPUT_DIR, 'stm-heldout-likelihood_{covariates}_{K}_{project}.tsv')
STM_EFFECT_TABLE=join(OUTPUT_DIR, 'stm-effect-table_{covariates}_{K}_{project}.tsv')
STM_EFFECT_PLOT=join(OUTPUT_DIR, 'stm-effect-table_{covariates}_{K}_{project}.png')
STM_HELDOUT_LIKELIHOOD_RATIO_FILE=join(OUTPUT_DIR, 'stm-heldout-ratio_{covariate_of_interest}_{covariates}_{K}_{project}.tsv')
STM_VALIDATE_LIKELIHOOD_RATIO_FILE=join(OUTPUT_DIR, 'stm-validate-ratio_{covariate_of_interest}_{covariates}_{K}_{project}.tsv')
COMBINED_HELDOUT_LIKELIHOOD_FILE=join(OUTPUT_DIR, 'stm-combined-heldout_{covariates}_{project}.tsv')
LIKELIHOOD_PLOT=join(OUTPUT_DIR,'likelihood-plot_{covariates1}_{covariates2}_{project}.pdf')
# Model Component Files
SIGMA_FILE=join(OUTPUT_DIR, 'stm-sigma_{covariates}_{K}_{project}.tsv')
GAMMA_FILE=join(OUTPUT_DIR, 'stm-gamma_{covariates}_{K}_{project}.tsv')

# Parameters
seed="123456"

rule convert_normalized_exposures_to_counts:
    input:
        STM_ALL_NORMALIZED_EXPOSURES_FILE,
        MC_FILE
    output:
        STM_ALL_COUNT_EXPOSURES_FILE
    script:
        "src/convert_normalized_to_count_exposures.py"

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
        seed
    input:
        TRAIN_MC_FILE,
        TEST_MC_FILE,
        TRAIN_FEATURE_FILE,
        TEST_FEATURE_FILE
    output:
        STM_TRAIN_NORMALIZED_EXPOSURES_FILE,
        STM_TEST_NORMALIZED_EXPOSURES_FILE,
        STM_TRAIN_SIGNATURES_FILE,
        STM_HELDOUT_LIKELIHOOD_RATIO_FILE,
    shell:
        '../src/tcsm_heldout_exposures.R "{input[0]}" "{input[1]}" "{input[2]}" "{input[3]}" {wildcards.K} {wildcards.covariates} {wildcards.covariate_of_interest} --seed={params[0]} --traine="{output[0]}" --teste="{output[1]}" --signature="{output[2]}" --heldout="{output[3]}"'

rule stm_heldout_exposures_validate:
    params:
        seed,
        FUNCTION_FILE
    input:
        MC_FILE,
        VALIDATE_MC_FILE,
        FEATURE_FILE,
        VALIDATE_FEATURE_FILE
    output:
        STM_NORMALIZED_EXPOSURES_FILE,
        STM_VALIDATE_NORMALIZED_EXPOSURES_FILE,
        STM_SIGNATURES_FILE,
        STM_VALIDATE_LIKELIHOOD_RATIO_FILE,
    shell:
        '../src/tcsm_heldout_exposures.R "{input[0]}" "{input[1]}" "{input[2]}" "{input[3]}" {wildcards.K} {wildcards.covariates} {wildcards.covariate_of_interest} --seed={params[0]} --traine="{output[0]}" --teste="{output[1]}" --signature="{output[2]}" --heldout="{output[3]}"'

rule run_TCSM_without_covariates:
    wildcard_constraints:
        covariates="NULL"
    params:
        seed
    input:
        MC_FILE
    output:
        STM_ALL_NORMALIZED_EXPOSURES_FILE,
        STM_EXOME_SIGNATURES_FILE,
        SIGMA_FILE,
    shell:
        'src/run_tcsm.R "{input[0]}" {wildcards.K} -s {params} --exposures="{output[0]}" --signatures="{output[1]}" --sigma="{output[2]}"'


rule run_TCSM_with_covariates:
    wildcard_constraints:
        covariates="(?!NULL).+"
    params:
        seed
    input:
        MC_FILE,
        FEATURE_FILE
    output:
        STM_ALL_NORMALIZED_EXPOSURES_FILE,
        STM_EXOME_SIGNATURES_FILE,
        STM_EFFECT_TABLE,
        SIGMA_FILE,
        GAMMA_FILE
    shell:
        'src/run_tcsm.R "{input[0]}" {wildcards.K} -c="{input[1]}" -s {params} --covariates {wildcards.covariates} --exposures="{output[0]}" --signatures="{output[1]}" --effect="{output[2]}" --sigma="{output[3]}" --gamma="{output[4]}"'

# calculate the likelihood of test samples after learning the model on train samples
rule tcsm_heldout_likelihood:
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
    shell:
        '../src/tcsm_heldout_likelihood.R "{input[0]}" "{input[1]}" {wildcards.K} --covariates={wildcards.covariates} --trainf="{input[2]}" --testf="{input[3]}" --heldout="{output[0]}"'
