# This is a assembly of functions which can pre-process data

# extract value data from bulk metadata
extract_value2list <- function(value) {
  ilist <- list()
  if (value == "Patient") {
    for (n in 2:ncol(bulk_clinical)) {
      a <- trimws(sapply(strsplit(as.vector(as.matrix(bulk_clinical[,n])), split = " ", fixed = T), '[' ,1))
      b <- trimws(sapply(strsplit(as.vector(as.matrix(bulk_clinical[,n])), split = " ", fixed = T), '[' ,3))
      ilist[[n-1]] <- b[which(a==value)]
      }
    } else { 
  for (n in 2:ncol(bulk_clinical)) {
    a <- trimws(sapply(strsplit(as.vector(as.matrix(bulk_clinical[2:nrow(bulk_clinical),n])), split = ":", fixed = T), '[' ,1))
    b <- trimws(sapply(strsplit(as.vector(as.matrix(bulk_clinical[2:nrow(bulk_clinical),n])), split = ":", fixed = T), '[' ,2))
    ilist[[n-1]] <- b[which(a==value)]
    }
  }
  return(flatten_chr(ilist))
}

# remove duplicate in single cell dataset
process_scdataset <- function(dataset) {
	a <- dataset %>% group_by(Symbol) %>% summarize(n= n()) %>% filter(n>1)
	idx <- list()
	for (n in 1:length(a$Symbol)) {
		b <- which(a$Symbol[n] == dataset$Symbol)[2]
		idx <- append(b, idx)
	}
	idx <- as.double(flatten_chr(idx))
	b <- dataset[!idx,]
  c <- as.matrix(b[,-c(1,2)])
	rownames(c) <- b$Symbol
   colnames(c) <- colnames(b)[-c(1,2)]

	return(c)
}

#  make median statistic
med <- function(data, indices) {
    resample <- data[indices]
    return(median(resample))
  }


zscore_transform <- function(x) {
  (x-mean(x))/sd(x)
}


make_genesignature <- function(dataset, log_dataset) {

	res2 <-  apply(log_dataset,2,mean)
	sumgsig <- summary(res2)
	names(res2)[which(res2<sumgsig[3])] -> low # low ID
	names(res2)[which(res2>=sumgsig[3])] -> high # high ID

	# Divide group
	dataset<- dataset%>% mutate(Novel_Signature = ifelse(ID %in% high, "High", "Low"))
	dataset$Novel_Signature <- factor(dataset$Novel_Signature)

	if (identical(dataset$ID, names(res2))) {
		dataset$Novel_Signature_score <- res2
	} else {
		res2 <- res2[dataset$ID]
		dataset$Novel_Signature_score <- res2
	}

	return(dataset)
}




