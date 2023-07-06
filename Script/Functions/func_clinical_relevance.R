
tpm <- function(counts, lengths) {
	rpk <- counts/(lengths/1000) # 1. rpk
	coef <- sum(rpk)/1e6 # 2. coef(per million scaling factor)
	rpk/coef # 3. tpm
}

zscore_transform <- function(x) {
  (x-mean(x))/sd(x)
}
