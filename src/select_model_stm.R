library(stm)
library(tm)


run.stm <- function(mutation.count.file, feature.file, covariates, seed, plot.output.file){
  mc.data <- read.delim(mutation.count.file, sep = '\t', header = TRUE, check.names=FALSE, row.names=1)
  mat <- data.matrix(mc.data)
  prep <- readCorpus(mat, type="dtm")
  feature.data <- read.delim(feature.file, sep = '\t', header = TRUE, row.names=1)

  pdf(plot.output.file)
  model.selection <- searchK(documents = prep$documents, vocab = prep$vocab,
     K=seq(2, 14), prevalence=covariates, heldout.seed=seed, data=feature.data, proportion=.9)
  print(model.selection)
  plot.searchK(model.selection)
  dev.off()
}


# this only works for CTM
mutation.count.file <- snakemake@input[[1]]
feature.file <- snakemake@input[[2]]
seed <- snakemake@params[[1]]
covariate.formula <- snakemake@params[[2]]
covariate.formula <- as.formula(covariate.formula)
plot.output.file <- snakemake@output[[1]]
run.stm(mutation.count.file, feature.file, covariate.formula, seed, plot.output.file)
