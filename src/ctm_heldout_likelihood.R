library(stm)
source(snakemake@params[[2]])
#source("stm_make_heldout.R")

run.ctm <- function(train.mutation.count.file, test.mutation.count.file, K, seed, heldout.performance.file){
  heldout <- make.heldout.obj(train.mutation.count.file, test.mutation.count.file)
  # the heldout object
  stm1 <- stm(documents=heldout$documents, vocab=heldout$vocab, K=K, seed=seed,
              max.em.its = 500, init.type = "Spectral")
  heldout.performance <- eval.heldout(stm1, heldout$missing)
  heldout.likelihood <- heldout.performance$expected.heldout
  df <- data.frame("likelihood" = heldout.likelihood, K)
  write.table(df, heldout.performance.file, sep="\t")
}

train.mutation.count.file <- snakemake@input[[1]]
test.mutation.count.file <- snakemake@input[[2]]
seed <- strtoi(snakemake@params[[1]])
K <- strtoi(snakemake@wildcards[["K"]])
heldout.performance.file <- snakemake@output[[1]]
run.ctm(train.mutation.count.file, test.mutation.count.file, K, seed, heldout.performance.file)
