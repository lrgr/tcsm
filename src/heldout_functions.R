library(stm)

load.stm.documents <- function(file.name){
  mc.data <- read.delim(file.name, sep = '\t', header = TRUE, check.names=FALSE, row.names=1)
  mat <- data.matrix(mc.data)
  corpus <- readCorpus(mat, type="dtm")
  prepDocuments(corpus$documents, corpus$vocab, lower.thresh=0)
}

# heldout is just a named list
#   $missing
#     $index
#     $docs - list where each element is a list containing the mutations for a given sample
#   $documents - list where each element is a list containing tthe mutations for a given sample
#   $vocab -
make.heldout.zero <- function(train.mutation.count.file, test.mutation.count.file){
  train.prep <- load.stm.documents(train.mutation.count.file)
  test.prep <- load.stm.documents(test.mutation.count.file)
  heldout <- list()
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

# Adapted from make.heldout function (see https://github.com/bstewart/stm/blob/410508f410a3bd53677dbe756f5afe8b424f13eb/R/heldout.R#L34)
make.heldout.obj <- function(train.mutation.count.file, test.mutation.count.file,
                         proportion=.5, seed=NULL) {
  if (proportion == 0){
    return(make.heldout.zero(train.mutation.count.file, test.mutation.count.file))
  }

  train.mc.data <- read.delim(train.mutation.count.file, sep = '\t', header = TRUE, check.names=FALSE, row.names=1)
  train.mat <- data.matrix(train.mc.data)
  test.mc.data <- read.delim(test.mutation.count.file, sep = '\t', header = TRUE, check.names=FALSE, row.names=1)
  test.mat <- data.matrix(test.mc.data)
  mat <- rbind(train.mat, test.mat)
  corpus <- readCorpus(mat, type="dtm")
  if(!is.null(seed)) set.seed(seed)

  # Convert the corpus to the internal STM format
  args <- asSTMCorpus(corpus$documents, corpus$vocab)
  documents <- args$documents
  vocab <- args$vocab

  index <- seq(nrow(train.mat)+1, nrow(train.mat)+nrow(test.mat))
  pie <- proportion
  missing <- vector(mode="list", length=length(index))
  ct <- 0
  for(i in index) {
    ct <- ct + 1
    doc <- documents[[i]]
    if(ncol(doc)<2) next
    doc <- rep(doc[1,], doc[2,])
    #how many tokens to sample? The max ensures at least one is sampled
    nsamp <- max(1,floor(pie*length(doc)))
    ho.index <- sample(1:length(doc), nsamp)
    tab <- tabulate(doc[ho.index])
    missing[[ct]] <- rbind(which(tab>0), tab[tab>0])
    tab <- tabulate(doc[-ho.index])
    documents[[i]] <- rbind(which(tab>0), tab[tab>0])
  }
  missing <- list(index=index, docs=missing)
  print(missing)
  #check the vocab
  indices <- sort(unique(unlist(lapply(documents, function(x) x[1,]))))

  #all sorts of nonsense ensues if there is missingness
  #first condition checks the vocab, second checks the documents
  if(length(indices)!=length(vocab) | any(unlist(lapply(missing$docs, is.null)))) {
    remove <- which(!(1:length(vocab)%in% indices))
    newind <- rep(0, length(vocab))
    newind[indices] <- 1:length(indices)
    new.map <- cbind(1:length(vocab), newind)
    #renumber the missing elements and remove 0's
    missing$docs <- lapply(missing$docs, function(d) {
      d[1,] <- new.map[match(d[1,], new.map[,1]),2]
      return(d[,d[1,]!=0, drop=FALSE])
    })
    #apply the same process to the documents
    documents <- lapply(documents, function(d) {
      d[1,] <- new.map[match(d[1,], new.map[,1]),2]
      return(d[,d[1,]!=0, drop=FALSE])
    })

    lens <- unlist(lapply(missing$docs, length))
    if(any(lens==0)) {
      missing$docs <- missing$docs[lens!=0]
      missing$index <- missing$index[lens!=0]
    }
    vocab <- vocab[indices]
  }
  #hooray.  return some stuff.
  heldout <- list(documents=documents,vocab=vocab, missing=missing)

  #you can get cases where these come out as non-integers...
  #recast everything just to be sure.
  heldout$documents <- lapply(heldout$documents, function(x) matrix(as.integer(x), nrow(x), ncol(x)))
  heldout$missing$docs <- lapply(heldout$missing$docs, function(x) matrix(as.integer(x), nrow(x), ncol(x)))
  class(heldout) <- "heldout"
  return(heldout)
}


get.heldout.ratio <- function(train.feature.file, test.feature.file, heldout, K, seed, covariates, covariate_of_interest){
  train.feature.data <- read.delim(train.feature.file, sep = '\t', header = TRUE, row.names=1)
  test.feature.data <- read.delim(test.feature.file, sep = '\t', header = TRUE, row.names=1)
  # print(train.feature.data)
  covariate.options <- sort(unique(train.feature.data[,covariate_of_interest]), decreasing=TRUE)
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
