#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library(stm))
source("../src/heldout_functions.R")


run.stm <- function(train.mutation.count.file, test.mutation.count.file, train.feature.file, test.feature.file, covariates, K, seed, heldout.performance.file){
  heldout <- make.heldout.obj(train.mutation.count.file, test.mutation.count.file)
  if (covariates == "NULL"){
    stm1 <- stm(documents=heldout$documents, vocab=heldout$vocab, K=K, seed=seed,
                max.em.its = 500, init.type = "Spectral")
  } else {
    train.feature.data <- read.delim(train.feature.file, sep = '\t', header = TRUE, row.names=1)
    test.feature.data <- read.delim(test.feature.file, sep = '\t', header = TRUE, row.names=1)
    feature.data <- rbind(train.feature.data, test.feature.data)
    covariate.formula <- as.formula(paste0("~", covariates))
    # the heldout object
    stm1 <- stm(documents=heldout$documents, vocab=heldout$vocab, K=K, seed=seed,
                prevalence = covariate.formula, max.em.its = 500, data=feature.data,
                init.type = "Spectral")
  }

  heldout.performance <- eval.heldout(stm1, heldout$missing)
  heldout.likelihood <- heldout.performance$expected.heldout
  df <- data.frame("likelihood" = heldout.likelihood, K)
  write.table(df, heldout.performance.file, sep="\t")
}

# create parser object
parser <- ArgumentParser()
# specify our desired options
parser$add_argument("trainmc", help="mutation count input file for training")
parser$add_argument("testmc", help="mutation count input file for test")
parser$add_argument("k", help="number of signatures to use")
parser$add_argument("--seed", help="random seed", default=123456)
parser$add_argument("--covariates", help="covariates (separated by +)", default="NULL")
parser$add_argument("--trainf", help="input file for covariate values for training samples")
parser$add_argument("--testf", help="input file for covariate values for test samples")
parser$add_argument("--heldout", default="heldout-log-likelihood.tsv", help="output file for heldout log-likelihood of test set")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

train.mutation.count.file <- args$trainmc
test.mutation.count.file <- args$testmc
train.feature.file <- args$trainf
test.feature.file <- args$testf
seed <- args$s
K <- strtoi(args$k)
covariates <- args$covariates
heldout.performance.file <- args$heldout


run.stm(train.mutation.count.file, test.mutation.count.file, train.feature.file, test.feature.file, covariates, K, seed, heldout.performance.file)
