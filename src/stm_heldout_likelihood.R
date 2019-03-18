library(stm)
source(snakemake@params[[2]])


run.stm <- function(train.mutation.count.file, test.mutation.count.file, train.feature.file, test.feature.file, covariates, K, seed, heldout.performance.file){
  heldout <- make.heldout.obj(train.mutation.count.file, test.mutation.count.file)
  train.feature.data <- read.delim(train.feature.file, sep = '\t', header = TRUE, row.names=1)
  test.feature.data <- read.delim(test.feature.file, sep = '\t', header = TRUE, row.names=1)
  feature.data <- rbind(train.feature.data, test.feature.data)
  # the heldout object
  stm1 <- stm(documents=heldout$documents, vocab=heldout$vocab, K=K, seed=seed,
              prevalence = covariates, max.em.its = 500, data=feature.data,
              init.type = "Spectral")
  heldout.performance <- eval.heldout(stm1, heldout$missing)
  heldout.likelihood <- heldout.performance$expected.heldout
  df <- data.frame("likelihood" = heldout.likelihood, K)
  write.table(df, heldout.performance.file, sep="\t")
}

train.mutation.count.file <- snakemake@input[[1]]
test.mutation.count.file <- snakemake@input[[2]]
train.feature.file <- snakemake@input[[3]]
test.feature.file <- snakemake@input[[4]]
seed <- strtoi(snakemake@params[[1]])
covariates <- snakemake@wildcards[["covariates"]]

covariate.formula <- as.formula(paste0("~", covariates))
print(covariate.formula)
K <- strtoi(snakemake@wildcards[["K"]])
heldout.performance.file <- snakemake@output[[1]]
run.stm(train.mutation.count.file, test.mutation.count.file, train.feature.file, test.feature.file, covariate.formula, K, seed, heldout.performance.file)
