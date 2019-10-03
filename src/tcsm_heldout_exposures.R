#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library(stm))
source("../src/heldout_functions.R")


run.stm <- function(train.mutation.count.file, test.mutation.count.file, train.feature.file, test.feature.file, covariates, K, seed, train.exposures.file, test.exposures.file, signature.output.file, covariate_of_interest, heldout.ratio.file){
  # run the model on just the training data
  train.prep <- load.stm.documents(train.mutation.count.file)
  train.feature.data <- read.delim(train.feature.file, sep = '\t', header = TRUE, row.names=1)
  covariate.formula <- as.formula(paste0("~", covariates))
  stm1 <- stm(documents=train.prep$documents, vocab=train.prep$vocab, K=K, seed=seed,
              prevalence = covariate.formula, max.em.its = 500, data=train.feature.data,
              init.type = "Spectral")
  # use fitNewDocuments with the test data and the predicted feature.data
  test.prep <- load.stm.documents(test.mutation.count.file)
  test.feature.data <- read.delim(test.feature.file, sep = '\t', header = TRUE, row.names=1)
  heldout <- make.heldout.obj(train.mutation.count.file, test.mutation.count.file)
  if (covariates != "NULL"){
    heldout.ratio <- get.heldout.ratio(train.feature.file, test.feature.file, heldout, K, seed, covariates, covariate_of_interest)
    df <- data.frame("heldout.ratio"=heldout.ratio)
    predicted.feature.data <- as.integer(heldout.ratio > 0)
    test.feature.data[covariates] <- predicted.feature.data
  } else {
    df <- data.frame(matrix(NA, nrow=nrow(test.feature.data), ncol=1))
  }
  df <- cbind(test.feature.data, df)
  write.table(df, heldout.ratio.file, sep="\t")
  test.results <- fitNewDocuments(model=stm1, documents=test.prep$documents, newData=test.feature.data,
                  origData=train.feature.data, prevalence=covariate.formula)
  test.exposures <- as.data.frame(test.results$theta)
  # print(test.exposures)
  # get the training exposure results
  train.exposures <- make.dt(stm1)
  train.exposures <- subset(train.exposures, select = -c(docnum))
  # get the sample names for the training data set
  train.sample.names <- rownames(read.delim(train.mutation.count.file, sep = '\t', header = TRUE, check.names=FALSE, row.names=1))
  # get the sample names for the test data set
  test.sample.names <- rownames(read.delim(test.mutation.count.file, sep = '\t', header = TRUE, check.names=FALSE, row.names=1))
  # set rows to be sample names
  rownames(train.exposures) <- train.sample.names
  rownames(test.exposures) <- test.sample.names
  colnames(test.exposures) <- colnames(train.exposures)
  #print(predicted.feature.data)
  #df <- data.frame("heldout.ratio"=heldout.ratio, "true.brcaness"=test.feature.data$feature1)
  write.table(train.exposures, train.exposures.file, sep="\t")
  write.table(test.exposures, test.exposures.file, sep="\t")
  mat <- stm1$beta$logbeta[[1]]
  signatures <- apply(mat, 1:2, exp)
  colnames(signatures) <- toupper(stm1$vocab)
  dt <- make.dt(stm1)
  rownames(signatures) <- colnames(dt)[-1]
  # save the signatures
  write.table(signatures, file=signature.output.file, sep="\t")
}

# create parser object
parser <- ArgumentParser()
# specify our desired options
parser$add_argument("trainmc", help="mutation count input file for training")
parser$add_argument("testmc", help="mutation count input file for test")
parser$add_argument("trainf", help="feature input file for training")
parser$add_argument("testf", help="feature input file for test")
parser$add_argument("k", help="number of signatures to use")
parser$add_argument("c", help="covariates (separated by +)")
parser$add_argument("coi", help="covariate of interest to be hidden in test")
parser$add_argument("--traine", help="train exposure output file")
parser$add_argument("--teste", help="test exposure output file")
parser$add_argument("--signature", help="output signature file")
parser$add_argument("--heldout", help="heldout output file")
parser$add_argument("--seed", help="random seed", default="123456")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

train.mutation.count.file <- args$trainmc
test.mutation.count.file <- args$testmc
train.feature.file <- args$trainf
test.feature.file <- args$testf
seed <- args$seed
K <- strtoi(args$k)
covariates <- args$c
covariate_of_interest <- args$coi
train.exposures.file <- args$traine
test.exposures.file <- args$teste
signatures.file <- args$signature
heldout.ratio.file <- args$heldout


run.stm(train.mutation.count.file, test.mutation.count.file, train.feature.file,
        test.feature.file, covariates, K, seed, train.exposures.file,
        test.exposures.file, signatures.file, covariate_of_interest, heldout.ratio.file)
