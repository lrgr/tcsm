library(stm)
source(snakemake@params[[2]])


run.stm <- function(train.mutation.count.file, test.mutation.count.file, train.feature.file, test.feature.file, covariates, K, seed, heldout.performance.file){
  heldout <- make.heldout.obj(train.mutation.count.file, test.mutation.count.file)
  heldout.ratio <- get.heldout.ratio(train.feature.file, test.feature.file, heldout, K, seed, covariates)
  test.feature.data <- read.delim(test.feature.file, sep = '\t', header = TRUE, row.names=1)
  df <- data.frame("heldout.ratio"=heldout.ratio)
  df <- cbind(test.feature.data, df)
  write.table(df, heldout.ratio.file, sep="\t")
}



train.mutation.count.file <- snakemake@input[[1]]
test.mutation.count.file <- snakemake@input[[2]]
train.feature.file <- snakemake@input[[3]]
test.feature.file <- snakemake@input[[4]]
seed <- strtoi(snakemake@params[[1]])
covariate.formula <- as.formula(paste0("~", snakemake@wildcards[["covariates"]]))
K <- strtoi(snakemake@wildcards[["K"]])
heldout.ratio.file <- snakemake@output[[1]]


run.stm(train.mutation.count.file, test.mutation.count.file, train.feature.file, test.feature.file, snakemake@wildcards[["covariates"]], K, seed, heldout.performance.file)
