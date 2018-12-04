library(stm)

# expected format is one column

run.ctm <- function(mutation.count.file, K, seed, exposure.output.file, signature.output.file){
  mc.data <- read.delim(mutation.count.file, sep = '\t', header = TRUE, check.names=FALSE, row.names=1)
  mat <- data.matrix(mc.data)
  corpus <- readCorpus(mat, type="dtm")
  prep <- prepDocuments(corpus$documents, corpus$vocab)
  stm1 <- stm(documents=prep$documents, vocab=prep$vocab, K=K, seed=seed,
              max.em.its = 500, init.type = "Spectral")

  # print(stm1$sigma)
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
  dt$docnum <- rownames(mc.data)
  # save the exposures for the training data set
  write.table(dt, file=exposure.output.file, sep="\t", row.names=FALSE)
}


# this only works for CTM
mutation.count.file <- snakemake@input[[1]]
seed <- snakemake@params[[1]]
K <- strtoi(snakemake@wildcards[["K"]])
exposure.output.file <- snakemake@output[[1]]
signature.output.file <- snakemake@output[[2]]
run.ctm(mutation.count.file, K, seed, exposure.output.file, signature.output.file)
