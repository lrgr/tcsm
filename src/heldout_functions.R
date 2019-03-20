library(stm)

load.stm.documents <- function(file.name){
  mc.data <- read.delim(file.name, sep = '\t', header = TRUE, check.names=FALSE, row.names=1)
  mat <- data.matrix(mc.data)
  corpus <- readCorpus(mat, type="dtm")
  prepDocuments(corpus$documents, corpus$vocab, lower.thresh=0)
}

# heldout is just a named list
#   $missing
#     $index
#     $docs - list where each element is a list containing the mutations for a given sample
#   $documents - list where each element is a list containing tthe mutations for a given sample
#   $vocab -
make.heldout.obj <- function(train.mutation.count.file, test.mutation.count.file){
  # need to pre-process test and train data sets together to avoid edge case where mutation categories don't line up
  train.mc.data <- read.delim(train.mutation.count.file, sep = '\t', header = TRUE, check.names=FALSE, row.names=1)
  test.mc.data <- read.delim(test.mutation.count.file, sep = '\t', header = TRUE, check.names=FALSE, row.names=1)
  mc.data <- rbind(train.mc.data, test.mc.data)
  mat <- data.matrix(mc.data)
  corpus <- readCorpus(mat, type="dtm")
  prep <- prepDocuments(corpus$documents, corpus$vocab, lower.thresh=0)
  n_training_samples = dim(train.mc.data)[1]
  n_test_samples = dim(test.mc.data)[1]
  heldout <- list()
  heldout$documents <- prep$documents[1:n_training_samples]
  heldout$missing <- list()
  heldout$missing$docs <- prep$documents[(n_training_samples+1):(n_training_samples+n_test_samples)]
  heldout$missing$index <- seq(n_training_samples+1, n_training_samples+n_test_samples)
  for (i in heldout$missing$index){
    heldout$documents[[i]] <- matrix(, nrow=2, ncol=0)
  }
  heldout$vocab <- prep$vocab
  heldout
}


get.heldout.ratio <- function(train.feature.file, test.feature.file, heldout, K, seed, covariates, covariate_of_interest){
  train.feature.data <- read.delim(train.feature.file, sep = '\t', header = TRUE, row.names=1)
  test.feature.data <- read.delim(test.feature.file, sep = '\t', header = TRUE, row.names=1)
  # print(train.feature.data)
  covariate.options <- sort(unique(train.feature.data[,covariate_of_interest]), decreasing=TRUE)
  stopifnot(length(covariate.options) == 2)
  test.brcad.data <- test.feature.data
  test.brcad.data[,covariate_of_interest] <- covariate.options[1]
  test.brcap.data <- test.feature.data
  test.brcap.data[,covariate_of_interest] <- covariate.options[2]
  brcad.feature.data <- rbind(train.feature.data, test.brcad.data)
  brcap.feature.data <- rbind(train.feature.data, test.brcap.data)
  covariate.formula <- as.formula(paste0("~", covariates))
  # the heldout object
  stm1 <- stm(documents=heldout$documents, vocab=heldout$vocab, K=K, seed=seed,
              prevalence = covariate.formula, max.em.its = 500, data=brcad.feature.data,
              init.type = "Spectral")
  brcad.heldout.performance <- eval.heldout(stm1, heldout$missing)
  stm1 <- stm(documents=heldout$documents, vocab=heldout$vocab, K=K, seed=seed,
              prevalence = covariate.formula, max.em.its = 500, data=brcap.feature.data,
              init.type = "Spectral")
  brcap.heldout.performance <- eval.heldout(stm1, heldout$missing)
  heldout.ratio <- brcad.heldout.performance$doc.heldout - brcap.heldout.performance$doc.heldout
}
