library(stm)

load.stm.documents <- function(file.name){
  mc.data <- read.delim(file.name, sep = '\t', header = TRUE, check.names=FALSE, row.names=1)
  mat <- data.matrix(mc.data)
  readCorpus(mat, type="dtm")
}

# heldout is just a named list
#   $missing
#     $index
#     $docs - list where each element is a list containing the mutations for a given sample
#   $documents - list where each element is a list containing tthe mutations for a given sample
#   $vocab -
make.heldout.obj <- function(train.mutation.count.file, test.mutation.count.file){
  train.prep <- load.stm.documents(train.mutation.count.file)
  test.prep <- load.stm.documents(test.mutation.count.file)
  heldout <- list("vocab")
  heldout$documents <- train.prep$documents
  heldout$missing <- list()
  heldout$missing$docs <- test.prep$documents
  heldout$missing$index <- seq(length(train.prep$documents)+1, length(train.prep$documents)+length(test.prep$documents))
  for (i in heldout$missing$index){
    heldout$documents[[i]] <- matrix(, nrow=2, ncol=0)
  }
  heldout$vocab <- train.prep$vocab
  heldout
}

get.heldout.ratio <- function(train.feature.file, test.feature.file, heldout, K, seed, covariates, covariate_of_interest){
  train.feature.data <- read.delim(train.feature.file, sep = '\t', header = TRUE, row.names=1)
  test.feature.data <- read.delim(test.feature.file, sep = '\t', header = TRUE, row.names=1)
  # print(train.feature.data)
  covariate.options <- sort(unique(train.feature.data[,covariate_of_interest]))
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
