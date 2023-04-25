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


 
 
