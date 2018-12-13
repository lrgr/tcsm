library(stm)


run.stm <- function(mutation.count.file, feature.file, covariates, K, seed, exposure.output.file, signature.output.file){
  mc.data <- read.delim(mutation.count.file, sep = '\t', header = TRUE, check.names=FALSE, row.names=1)
  mat <- data.matrix(mc.data)
  corpus <- readCorpus(mat, type="dtm")
  prep <- prepDocuments(corpus$documents, corpus$vocab)
  feature.data <- read.delim(feature.file, sep = '\t', header = TRUE, row.names=1)
  covariate.formula <- as.formula(paste0("~", covariates))
  stm1 <- stm(documents=prep$documents, vocab=prep$vocab, K=K, seed=seed, prevalence = covariate.formula, max.em.its = 500, data=feature.data, init.type = "Spectral")
  effect <- estimateEffect(covariate.formula, stm1, metadata=feature.data)
  effect.summary <- summary(effect)
  print(effect.summary)
  effect.tables <- effect.summary$tables
  results <- lapply(effect.summary$tables, function(x) x[, "Estimate"])
  effect.frame <- as.data.frame(do.call(rbind, results))
  write.table(effect.frame, file=snakemake@output[[3]], sep="\t")
  # print(stm1$mu)
  # print(stm1$sigma)
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


# this only works for STM
mutation.count.file <- snakemake@input[[1]]
feature.file <- snakemake@input[[2]]
seed <- snakemake@params[[1]]
K <- strtoi(snakemake@wildcards[["K"]])
exposure.output.file <- snakemake@output[[1]]
signature.output.file <- snakemake@output[[2]]
run.stm(mutation.count.file, feature.file, snakemake@wildcards[["covariates"]], K, seed, exposure.output.file, signature.output.file)
#
# if (nchar(covariate.formula) == 0){
#   mu.output.file <- snakemake@output[[4]]
#   sigma.output.file <- snakemake@output[[5]]
#   print(K)
#   run.ctm.cv(train.input.file, test.input.file, signature.output.file, train.exposure.output.file, test.exposure.output.file, mu.output.file, sigma.output.file, K)
# } else {
#   covariate.formula <- as.formula(covariate.formula)
#   run.stm.cv(train.input.file, test.input.file, covariate.formula, signature.output.file, train.exposure.output.file, test.exposure.output.file, K)
# }
