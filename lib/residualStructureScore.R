positiveNegativeApproxError <- function(m, c=2){
  positive_m <- m
  positive_m[positive_m<=0] <- 0.00000001 #avoid all zero rows
  nmf_result <- nmf(positive_m,c,method='lee')
  approx_positive_m <- (nmf_result@fit@W) %*% (nmf_result@fit@H)
  
  negative_m <- -m
  if(sum(negative_m<=0)==ncol(m)^2){
    approx_negative_m <- 0
  }else{
    negative_m[negative_m<0] <- 0.00000001
    nmf_result <- nmf(negative_m,c,method='lee')
    approx_negative_m <- (nmf_result@fit@W) %*% (nmf_result@fit@H)
  }
  
  return(c(1-sqrt(sum((approx_positive_m-positive_m)^2))/sqrt(sum(positive_m^2)),
           1-sqrt(sum((approx_negative_m-negative_m)^2))/sqrt(sum(negative_m^2))))
}