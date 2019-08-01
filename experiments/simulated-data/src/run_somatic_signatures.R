suppressPackageStartupMessages(library("SomaticSignatures"))

set.seed(strtoi(snakemake@params[["seed"]]))
df <- read.delim(snakemake@input[[1]], sep = '\t', header = TRUE, check.names=FALSE)
# df is samples-by-mutation-categories
# SomaticSignatures expects categories-by-samples
mat <- t(data.matrix(df))
mat <- mat[-1,]
colnames(mat) <- seq(0, ncol(mat)-1)
mat <- mat[as.logical(rowSums(mat != 0)), ]
out <- identifySignatures(mat, strtoi(snakemake@wildcards[["K"]]))
# their signatures matrix is categories-by-signatures
# our expected format for signatures is signatures-by-categories
signatures <- t(signatures(out))
write.table(signatures, file=snakemake@output[[1]], sep="\t")
# their (normalized) exposures matrix is samples-by-signatures (same as our expected format)
exposures <- samples(out)
# normalize the exposures so samples sum to 1
normalized_exposures <- apply(exposures, 1, function(i) i/sum(i))
# exposures <- scale(exposures, scale=rowSums(exposures))
normalized_exposures <- t(normalized_exposures)
write.table(normalized_exposures, file=snakemake@output[[2]], sep="\t")
