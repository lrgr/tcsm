#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library(stm))


run.tcsm <- function(mutation.count.file, feature.file, covariates, K, seed, exposure.output.file, signature.output.file, effect.output.file, sigma.output.file, gamma.output.file){
  #covariate.list <- strsplit(covariates, "\\+")[[1]]
  #print(c("default", covariate.list))
  mc.data <- read.delim(mutation.count.file, sep = '\t', header = TRUE, check.names=FALSE, row.names=1)
  mat <- data.matrix(mc.data)
  corpus <- readCorpus(mat, type="dtm")
  prep <- prepDocuments(corpus$documents, corpus$vocab)
  feature.data <- read.delim(feature.file, sep = '\t', header = TRUE, row.names=1)
  covariate.formula <- as.formula(paste0("~", covariates))
  stm1 <- stm(documents=prep$documents, vocab=prep$vocab, K=K, seed=seed,
              prevalence = covariate.formula, max.em.its = 500, data=feature.data,
              init.type = "Spectral", sigma=0)
  effect <- estimateEffect(covariate.formula, stm1, metadata=feature.data)
  effect.summary <- summary(effect)
  print(effect.summary)
  tc <- topicCorr(stm1, method="simple")
  # print(tc)
  # plot.topicCorr(tc)
  effect.tables <- effect.summary$tables
  results <- lapply(effect.summary$tables, function(x) x[, "Estimate"])
  effect.frame <- as.data.frame(do.call(rbind, results))
  write.table(effect.frame, file=effect.output.file, sep="\t")
  covariate.list <- strsplit(covariates, "\\+")[[1]]
  gamma <- stm1$mu$gamma
  print(gamma)
  if (covariates != "NULL"){
    rownames(gamma) <- c("default", covariate.list)
  }
  write.table(gamma, file=gamma.output.file, sep="\t")
  write.table(stm1$sigma, file=sigma.output.file, sep="\t")
  # plot.estimateEffect(prep)
  # searchK(documents = train.out$documents, vocab = train.out$vocab, K=seq(3, 10))
  # process the signatures and save them
  # the K-by-V matrix logbeta contains the natural log of the probability of seeing each word conditional on the topic
  mat <- stm1$beta$logbeta[[1]]
  signatures <- apply(mat, 1:2, exp)
  colnames(signatures) <- toupper(stm1$vocab)
  dt <- make.dt(stm1)
  rownames(signatures) <- colnames(dt)[-1]
  # save the signatures
  write.table(signatures, file=signature.output.file, sep="\t")
  # save the exposures for the training data set
  dt$docnum <- rownames(mc.data)
  # print(colnames(dt)[-1])
  # colnames(dt)[-1] <- rownames(signatures)
  write.table(dt, file=exposure.output.file, sep="\t", row.names=FALSE)
}


# create parser object
parser <- ArgumentParser()
# specify our desired options
parser$add_argument("-m", help="mutation count input file")
parser$add_argument("-c", help="covariate input file")
parser$add_argument("-e", "--exposures", help="normalized exposure output file")
parser$add_argument("--signatures", help="exome signature output file")
parser$add_argument("--effect", help="effect output file")
parser$add_argument("--sigma", help="sigma output file")
parser$add_argument("--gamma", help="gamma output file")
parser$add_argument("-s", help="random seed")
parser$add_argument("--covariates", help="covariates (separated by +)")
parser$add_argument("-k", help="number of signatures to use")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

mutation.count.file <- args$m
feature.file <- args$c
seed <- args$s
K <- strtoi(args$k)
exposure.output.file <- args$exposures
signature.output.file <- args$signatures
effect.output.file <- args$effect
sigma.output.file <- args$sigma
gamma.output.file <- args$gamma
covariates <- args$covariates


# mutation.count.file <- snakemake@input[[1]]
# feature.file <- snakemake@input[[2]]
#
# seed <- snakemake@params[[1]]
# K <- strtoi(snakemake@wildcards[["K"]])
# exposure.output.file <- snakemake@output[[1]]
# signature.output.file <- snakemake@output[[2]]
# effect.output.file <- snakemake@output[[3]]
# sigma.output.file <- snakemake@output[[4]]
# gamma.output.file <- snakemake@output[[5]]
run.tcsm(mutation.count.file, feature.file, covariates, K, seed, exposure.output.file, signature.output.file, effect.output.file, sigma.output.file, gamma.output.file)
