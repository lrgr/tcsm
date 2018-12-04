library(stm)
source(snakemake@params[[3]])
#source("stm_make_heldout.R")

run.stm <- function(train.mutation.count.file, test.mutation.count.file, train.feature.file, test.feature.file, covariates, K, seed, train.exposures.file, test.exposures.file){
  heldout <- make.heldout.obj(train.mutation.count.file, test.mutation.count.file)
  heldout.ratio <- get.heldout.ratio(train.feature.file, test.feature.file, heldout, K, seed, covariates)
  predicted.feature.data <- as.integer(heldout.ratio > 0)
  # run the model on just the training data
  train.prep <- load.stm.documents(train.mutation.count.file)
  train.feature.data <- read.delim(train.feature.file, sep = '\t', header = TRUE, row.names=1)

  stm1 <- stm(documents=train.prep$documents, vocab=train.prep$vocab, K=K, seed=seed,
              prevalence = covariates, max.em.its = 500, data=train.feature.data,
              init.type = "Spectral")
  # use fitNewDocuments with the test data and the predicted feature.data
  test.prep <- load.stm.documents(test.mutation.count.file)
  #print(test.prep$documents)
  #print(predicted.feature.data)
  #print(train.feature.data)
  test.feature.data <- read.delim(test.feature.file, sep = '\t', header = TRUE, row.names=1)
  test.feature.data$feature1 <- predicted.feature.data
  test.results <- fitNewDocuments(model=stm1, documents=test.prep$documents, newData=test.feature.data,
                  origData=train.feature.data, prevalence=covariates)
  test.exposures <- as.data.frame(test.results$theta)
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
}



train.mutation.count.file <- snakemake@input[[1]]
test.mutation.count.file <- snakemake@input[[2]]
train.feature.file <- snakemake@input[[3]]
test.feature.file <- snakemake@input[[4]]
seed <- strtoi(snakemake@params[[1]])
covariate.formula <- as.formula(snakemake@params[[2]])
K <- strtoi(snakemake@wildcards[["K"]])
train.exposures.file <- snakemake@output[[1]]
test.exposures.file <- snakemake@output[[2]]


run.stm(train.mutation.count.file, test.mutation.count.file, train.feature.file,
        test.feature.file, covariate.formula, K, seed, train.exposures.file,
        test.exposures.file)
