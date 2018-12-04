library(stm)
# library(tm)

run.stm <- function(train.mutation.count.file, test.mutation.count.file, K, seed, heldout.performance.file){
  train.mc.data <- read.delim(train.mutation.count.file, sep = '\t', header = TRUE, check.names=FALSE, row.names=1)
  test.mc.data <- read.delim(test.mutation.count.file, sep = '\t', header = TRUE, check.names=FALSE, row.names=1)
  #print(mc.data)
  train.mat <- data.matrix(train.mc.data)
  test.mat <- data.matrix(test.mc.data)
  train.prep <- readCorpus(train.mat, type="dtm")
  test.prep <- readCorpus(test.mat, type="dtm")
  # the heldout object
  stm1 <- stm(documents=train.prep$documents, vocab=train.prep$vocab, K=K, seed=seed,
              max.em.its = 500, init.type = "Spectral")
  new <- fitNewDocuments(model=stm1, documents=test.prep$documents,
    returnPosterior=TRUE)
  # print(new)
  # we need to get the number of tokens per document
  # print(test.prep$documents)
  ntokens <- unlist(lapply(test.prep$documents, function(x) sum(x[2,])))
  # print(heldout.performance$ntokens)
  # print(new$bound/heldout.performance$ntokens)
  bound <- mean(new$bound/ntokens)

  df <- data.frame(bound, K)
  write.table(df, heldout.performance.file, sep="\t")
  # print(cor(new$bound/heldout.performance$ntokens, heldout.performance$doc.heldout))
  # #print(prep)
  # print(heldout.performance$expected.heldout)
  #print(heldout.performance$doc.heldout)
  # saveRDS(new, file=heldout.performance.file)
}


train.mutation.count.file <- snakemake@input[[1]]
test.mutation.count.file <- snakemake@input[[2]]
seed <- snakemake@params[[1]]
K <- strtoi(snakemake@wildcards[["K"]])
heldout.performance.file <- snakemake@output[[1]]
run.stm(train.mutation.count.file, test.mutation.count.file, K, seed, heldout.performance.file)
