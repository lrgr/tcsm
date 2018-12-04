library(stm)
library(tm)


run.ctm <- function(input.file, seed, plot.output.file){
  mc.data <- read.delim(mutation.count.file, sep = '\t', header = TRUE, check.names=FALSE, row.names=1)
  mat <- data.matrix(mc.data)
  prep <- readCorpus(mat, type="dtm")
  pdf(plot.output.file)
  model.selection <- searchK(documents = prep$documents, vocab = prep$vocab,
    K=seq(2, 14), heldout.seed=seed, proportion=.9)
  print(model.selection)
  plot.searchK(model.selection)
  dev.off()
}


# this only works for CTM
mutation.count.file <- snakemake@input[[1]]
seed <- snakemake@params[[1]]
plot.output.file <- snakemake@output[[1]]
run.ctm(mutation.count.file, seed, plot.output.file)
