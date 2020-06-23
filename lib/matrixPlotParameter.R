# ----- Define a function for plotting a matrix ----- #
matrixPlotParameter <- function(x, ...){
  min <- min(x)
  max <- max(x)
  #abs<-6.9
  #min<--abs
  #min<--0.78
  #max<-20
  #min<-0
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  par(oma = c(1,1,2,1))
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(5,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- colorRampPalette(rev(brewer.pal(8,"RdBu")))(length(xLabels)*length(yLabels))
  ColorRamp <- viridis::viridis(n = length(xLabels)*length(yLabels))
  # ColorRamp <- two.colors(n = length(xLabels)*length(yLabels),end = rgb(0,0,0),
  #                         start = rgb(255,255,255,maxColorValue = 255) )
  #cm.colors(256) #rgb( seq(0,1,length=256),  # Red
  #      seq(0,1,length=256),  # Green
  #     seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp)) #min max
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  #yLabels <- yLabels[reverse]
  #x <- x[reverse,]
  
  # Data Map
  par(mar = c(1,3,1,1),mgp=c(0,0.7,0))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=1.3, tck=-0.015)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=1.3, tck=-0.015)
  mtext(side=3, text=title, line=0.8, cex=1.6, font=2)
  # Color Scale
  par(mar = c(1,2,1,1))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n", cex.axis=1.2)
  
  layout(1)
}

# ----- Define a function for plotting a contour plot ----- #
my.image.plot <- function(..., add = FALSE, xlab=NULL, ylab=NULL, my.main=NULL, cex.main=1, cex.lab = 1,
                          xLabels=NULL, yLabels=NULL,
                          breaks= NULL, nlevel = 64, col = NULL,  
                          horizontal = FALSE, legend.shrink = 0.9, legend.width = 1.2, 
                          legend.mar = ifelse(horizontal, 3.1, 5.1), legend.lab = NULL,
                          legend.line= 2,                    
                          graphics.reset = FALSE, bigplot = NULL, smallplot = NULL, 
                          legend.only = FALSE,  lab.breaks = NULL, 
                          axis.args = NULL, legend.args = NULL, legend.cex=1.0, midpoint = FALSE, border = NA, 
                          lwd = 1, verbose=FALSE) {
  # Thanks to S. Koehler and  S. Woodhead
  # for comments on making this a better function
  #
  # save current graphics settings
  old.par <- par(no.readonly = TRUE)
  # set defaults for color scale 
  # note this works differently than the image function.
  if( is.null(col))  {
    #col<-  tim.colors(nlevel)}
    #col <- two.colors(nlevel,start = rgb(0,0,0),end=rgb(255,255,255,maxColorValue = 255))}
    col <- viridis::viridis(nlevel)}
  else{
    nlevel<- length( col)
  }
  #  figure out zlim from passed arguments
  #  also set the breaks for colors if they have not been passed, 
  info <- imagePlotInfo(..., breaks=breaks, nlevel=nlevel)
  # breaks have been computed if not passed in the call
  breaks<- info$breaks
  if( verbose){
    print(info)
  }
  if (add) {
    big.plot <- old.par$plt
  }
  if (legend.only) {
    graphics.reset <- TRUE
  }
  if (is.null(legend.mar)) {
    legend.mar <- ifelse(horizontal, 3.1, 5.1)
  }
  # figure out how to divide up the plotting real estate 
  temp <- imageplot.setup(add = add, legend.shrink = legend.shrink, 
                          legend.width = legend.width, legend.mar = legend.mar, 
                          horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
  # bigplot has plotting region coordinates for image
  # smallplot has plotting coordinates for legend strip
  smallplot <- temp$smallplot
  bigplot <- temp$bigplot
  # draw the image in bigplot, just call the R base function
  # or poly.image for polygonal cells
  # note the logical switch
  # for poly.grid is parsed out of call from image.plot.info
  par(oma=c(3.8,4,2.8,1))
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(5.5,1), heights=c(1,1))
  if (!legend.only) {
    if (!add) {
      #par(plt = bigplot)
      par(mar=rep(0,4), mgp=c(0,0.8,0))
    }
    if (!info$poly.grid) {
      image(..., breaks=breaks, add = add, col = col, ylab="", xlab="", axes=FALSE)
      mtext(side=2, text=ylab, line=2.2, las=2, cex=cex.lab)
      mtext(side=1, text=xlab, line=2.7, cex=cex.lab)
      mtext(side=3, text=my.main, line=1, cex=cex.main, font=2)
      if(!is.null(xLabels)){
        axis(BELOW<-1, at=c(0.05,0.95), labels=xLabels, cex.axis=1.8, tck=-0.015)
        axis(LEFT <-2, at=c(0.05,0.95), labels=yLabels, las= HORIZONTAL<-1, cex.axis=1.8, tck=-0.015)
      }
    }
    else {
      poly.image(..., add = add, col = col, midpoint = midpoint, 
                 border = border, lwd.poly = lwd)
    }
    big.par <- par(no.readonly = TRUE)
  }
  ##
  ## check dimensions of smallplot
  # if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
  #   par(old.par)
  #   stop("plot region too small to add legend\n")
  # }
  # Following code draws the legend using the image function
  # and a one column image.
  # What might be confusing is the values of the "image" are the same 
  # as the locations on the legend axis.
  # Moreover the image values are in the middle of each breakpoint category
  # thanks to Tobias Nanu Frechen and Matthew Flickinger 
  # for sorting out some problems with the breaks position in the legend.
  ix <- 1
  iy<- breaks
  nBreaks<- length( breaks)
  midpoints<- (breaks[1:(nBreaks-1)] +  breaks[2:nBreaks] )/2
  iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints)) 
  if( verbose){print(breaks)
    print( midpoints)
    print( ix)
    print( iy)
    print( iz)
    print( col)}  
  #
  # next par call sets up a new plotting region just for the legend strip
  # at the smallplot coordinates
  par(mar=c(0,1.5,0,2), mgp=c(0,0.45,0))
  # draw color scales the two  cases are horizontal/vertical 
  # add a label if this is passed.
  if (!horizontal) {
    image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
          ylab = "", col = col, breaks=breaks, axes=FALSE)
    axis(4, at=c(0,0.6,1.2,1.8,2.4), labels=c(0,0.6,1.2,1.8,2.4), las= HORIZONTAL<-1, 
         cex.axis=1.6, tck=-0.15)
  }
  else {
    image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
          ylab = "", col = col, breaks=breaks, axes=FALSE)
    axis(4, at=c(0,0.6,1.2), labels=c(0,0.6,1.2), las= HORIZONTAL<-1, 
         cex.axis=1.6, tck=-0.1)
  }
  
  box()
  #
  # add a label to the axis if information has been  supplied
  # using the mtext function. The arguments to mtext are
  # passed as a list like the drill for axis (see above)
  #
  if (!is.null(legend.lab)) {
    legend.args <- list(text = legend.lab, side = ifelse(horizontal, 
                                                         1, 4), 
                        line = legend.line, cex=legend.cex)
    #                    just guessing at a good default for line argument!
  }
  # add the label using mtext function
  if (!is.null(legend.args)) {
    do.call(mtext, legend.args)
  }
  #
  # clean up graphics device settings
  # reset to larger plot region with right user coordinates.
  mfg.save <- par()$mfg
  if (graphics.reset | add) {
    par(old.par)
    par(mfg = mfg.save, new = FALSE)
    invisible()
  }
  else {
    par(big.par)
    par(plt = big.par$plt, xpd = FALSE)
    par(mfg = mfg.save, new = FALSE)
    # Suggestion from Karline Soetaert <Karline.Soetaert@nioz.nl>
    # this is to reset margins to be based on the mar arguments
    #      par(mar = par("mar"))  or
    #      par(mar = big.par$mar)
    # unfortunately this causes problems by allowing plotting outside of the
    # original plot region.
    invisible()
  }
}

# ----- Define a function for plotting multiple matrices ----- #
myMultiMatrixPlot <- function(X,xLabels,yLabels,min,max,colorPalette,no_matrix,n_row,
                              if.emp=TRUE,emp.ROW,emp.COL,emp.color,title_lst, reverse.col=TRUE,
                              title_cex=2.8, axis_cex=2.2, col_axis, fake_matrix=TRUE, out.title=NULL, 
                              matrix.mar=NULL, legend.mar=NULL, title.line=NULL){
  if(is.null(matrix.mar)){
    if(is.null(out.title)){
      par(oma=c( 0, 1, 0, 0 )+0.5,mgp=c(8,1,0))
    }else{
      par(oma=c( 0.5, 1, 2.5, 0 )+0.5,mgp=c(0.5,1,0))
    }
  }else{
    par(oma=c( 0.5, 0.5, 3, 2 ),mgp=c(0,1,0))
  }
  
  if(no_matrix==2){
    layout(matrix(c(1:(no_matrix+2)), nrow=1), 
           widths=c(rep(c(6,1.7),no_matrix-1),6,1.2), heights=1)
  }else if(no_matrix==3){
    layout(matrix(c(1:(no_matrix+3)), nrow=1), 
           widths=c(rep(c(6,1.7),no_matrix-1),6,1.2), heights=1)
  }else{
    layout(cbind(matrix(data=c(1:no_matrix), nrow=n_row, ncol=no_matrix/n_row, byrow = TRUE),
                 rep(no_matrix+1,n_row)), 
           widths=c(rep(6,no_matrix/n_row),1), heights=rep(1,n_row))
  }
  for(i in 1:no_matrix){
    x <- X[[i]]
    
    if(is.list(yLabels)){
      yLabels_temp <- yLabels[[i]]
      xLabels_temp <- xLabels[[i]]
      min_temp <- min[[i]]
      max_temp <- max[[i]]
    }else{
      yLabels_temp <- yLabels
      xLabels_temp <- xLabels
      min_temp <- min
      max_temp <- max
    }
    
    if(is.list(colorPalette)){
      ColorRamp <- colorRampPalette(brewer.pal(8,colorPalette[[i]]))(length(xLabels_temp)*length(yLabels_temp))
      #ColorRamp <- colorRampPalette(brewer.pal(8,colorPalette[[i]]))(length(xLabels_temp)*length(yLabels_temp))
    }else{
      if(reverse.col){
        ColorRamp <- viridis::viridis(length(xLabels_temp)*length(yLabels_temp))
        ColorRamp <- colorRampPalette(rev(brewer.pal(8,colorPalette)))(length(xLabels_temp)*length(yLabels_temp))
      }else{
        ColorRamp <- viridis::viridis(length(xLabels_temp)*length(yLabels_temp))
        ColorRamp <- colorRampPalette(brewer.pal(8,colorPalette))(length(xLabels_temp)*length(yLabels_temp))
      }
    }
    ColorLevels <- seq(min_temp, max_temp, length=length(ColorRamp)) 
    
    # Reverse Y axis
    # reverse <- nrow(x) : 1
    # yLabels_temp <- yLabels_temp[reverse]
    # x <- x[reverse,]
    
    # Data Map
    if(is.null(matrix.mar)){
      par(mar = c(axis_cex,axis_cex,title_cex,1)*1.1)
    }else{
      par(mar = matrix.mar)
    }
    image(1:length(xLabels_temp), 1:length(yLabels_temp), t(x), col=ColorRamp, xlab="",
          ylab="", axes=FALSE, zlim=c(min_temp,max_temp))
    if(is.null(title.line)){
      mtext(3, text=title_lst[[i]], cex=title_cex)
      
      axis(BELOW<-1, at=1:length(xLabels_temp), labels=xLabels_temp, cex.axis=axis_cex)
      axis(LEFT <-2, at=1:length(yLabels_temp), labels=yLabels_temp, las= HORIZONTAL<-1, cex.axis=axis_cex)
    }else{
      mtext(3, text=title_lst[[i]], cex=title_cex, line=title.line$title, font=2)
      mtext(text=xLabels_temp, side=3, at=1:length(xLabels_temp), cex=axis_cex*0.65, line=title.line$top)
      #mgp.axis(3,mgp=c(0,title.line$top,0),at=1:length(xLabels_temp), labels=xLabels_temp,
      #         cex.axis=axis_cex, lwd.ticks=0.2, lwd=0)
      mgp.axis(2,mgp=c(0,title.line$left,0),at=1:length(yLabels_temp), labels=yLabels_temp, 
               las= HORIZONTAL<-1, cex.axis=axis_cex, lwd=0)
      #axis(BELOW<-1, at=1:length(xLabels_temp), labels=xLabels_temp, cex.axis=axis_cex, line=title.line$bottom)
      #axis(LEFT <-2, at=1:length(yLabels_temp), labels=yLabels_temp, las= HORIZONTAL<-1, cex.axis=axis_cex, line=title.line$left)
    }
    
    if(if.emp){
      #Emphasize on some entries
      poly <- matrix.poly(1:length(xLabels_temp),1:length(yLabels_temp),t(x),ROW=emp.COL[[i]],COL=length(yLabels_temp)-emp.ROW[[i]]+1)
      polygon(poly,col=ifelse(emp.color[[i]]=="transparent",rgb(1,1,1,0),emp.color[[i]]),border=1)
    }
    if(i<no_matrix&fake_matrix){
      plot(0,0,xlim=c(0,5),ylim=c(0,5),type="n", xlab="", ylab="", axes = FALSE, frame=FALSE)
    }
  }
  # Color Scale
  if(is.null(legend.mar)){
    par(mar = c(axis_cex,0.6,title_cex,axis_cex)) 
  }else{
    par(mar = legend.mar,mgp=c(0,0.8,0))
  }
  
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        axes=FALSE)
  axis(RIGHT <- 4, at=col_axis, las=2,
       labels=col_axis,
       cex.axis=axis_cex*0.9)
  mtext(out.title, cex=2.4, outer=TRUE, line=-0.1, at=0.45, font=2)
  layout(1)
}
