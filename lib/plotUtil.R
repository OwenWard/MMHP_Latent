add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

cap.boxplot <- function(values, quantile=c(0.05,0.95)){
  return.values <- values
  return.values[values>quantile(values,quantile[2])] <- quantile(values,quantile[2])
  return.values[values<quantile(values,quantile[1])] <- quantile(values,quantile[1])
  return(return.values)
}

makeTransparent <- function(..., alpha=0.5) {
  
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
  
}

cRamp <- function(x, if.reverse = FALSE, bound = NULL){
  if(is.null(bound)){
    if(!if.reverse){
      cols <- colorRamp(brewer.pal(9,"Blues"))((x-min(x))/diff(range(x)))
    }else{
      cols <- colorRamp(brewer.pal(9,"Blues"))((max(x)-x)/diff(range(x)))
    }
  }else{
    bound[1] <- min(min(x),bound[1]) 
    bound[2] <- max(max(x),bound[2]) 
    if(!if.reverse){
      cols <- colorRamp(brewer.pal(9,"Blues"))((x-bound[1])/diff(bound))
    } else{
      cols <- colorRamp(brewer.pal(9,"Blues"))((bound[2]-x)/diff(bound))
    }
  }
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  


cLegend <- function(x, if.reverse = FALSE){
  max.x <- max(x)
  min.x <- min(x)
  x.vec <- seq(min.x,max.x,length=100)
  if(if.reverse){
    cols <- colorRamp(rev(brewer.pal(9,"Blues")))(x.vec/(max.x-min.x))
  }else{
    cols <- colorRamp(brewer.pal(9,"Blues"))(x.vec/(max.x-min.x))
  }
  return(list(x.vec=x.vec, 
              col.vec=apply(cols, 1, function(xt) rgb(xt[1], xt[2], xt[3], maxColorValue=255))))
}  

matrixToVecIdx <- function(i,j,N,byCol=TRUE){
  if(byCol){
    return(N*(j-1)+i)
  }else{
    return(N*(i-1)+j)
  }
}