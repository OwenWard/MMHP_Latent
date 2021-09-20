addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=TRUE)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

myCircularPlot <- function(cur.matrix, cur.order, cur.color,
                           cur.gap, cur.ID = c(1:12)){
  par(mfrow = c(1, 1),
      mgp = c(0, 0, 0),
      mar = c(1, 0, 0, 0),
      oma = c(0, 0, 0, 0))
  df1 <- data.frame(order = cur.order,
                    ID = cur.ID,
                    stringsAsFactors = FALSE)
  if(length(cur.order)==12){
    m <- cur.matrix
    df1 <- data.frame(order = cur.order,
                      ID = as.character(c(1:12)),
                      stringsAsFactors = FALSE)
    
  }else{
    m <- cur.matrix
    df1 <- data.frame(order = cur.order,
                      ID = cur.ID,
                      stringsAsFactors = FALSE)
  }
  dimnames(m) <- list(orig = df1$ID, dest = df1$ID)
  ### Sort order of data.frame and matrix for plotting in circos
  df1 <- arrange(df1, order)
  df1$ID <- factor(df1$ID, levels = df1$ID)
  m <- m[levels(df1$ID), levels(df1$ID)]
  ### Define ranges of circos sectors and their colors (both of the sectors and the links)
  df1$xmin <- 0
  df1$xmax <- rowSums(m) + colSums(m)
  n <- nrow(df1)
  
  df1$rcol<- cur.color#rgb(df1$r, df1$g, df1$b, max = 255)
  df1$lcol<- addalpha(df1$rcol, alpha = rep(0.5, n)) #rgb(df1$r, df1$g, df1$b, alpha=100, max = 255)
  
  ### Plot sectors (outer part)
  circos.clear()
  
  ### Basic circos graphic parameters
  circos.par(
    cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0.15),
    start.degree = 90,
    gap.degree = cur.gap
  )

  ### Sector details
  circos.initialize(factors = df1$ID, xlim = cbind(df1$xmin, df1$xmax))
  ### Plot sectors
  circos.trackPlotRegion(
    ylim = c(0, 1), factors = df1$ID, track.height = 0.1,
    # panel.fun for each sector
    panel.fun = function(x, y) {
      # select details of current sector
      name <- get.cell.meta.data("sector.index")
      i <- get.cell.meta.data("sector.numeric.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")

      # text direction (dd) and adjusmtents (aa)
      theta <- circlize(mean(xlim), 1.3)[1, 1] %% 360
      dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
      aa <- c(1, 0.5)
      if (theta < 90 || theta > 270) aa <- c(0, 0.5)

      # plot ID labels
      circos.text(x = mean(xlim), y = 1.7, labels = name, facing = dd, cex = 2.6, adj = aa)

      # plot main sector
      circos.rect(
        xleft = xlim[1], ybottom = ylim[1], xright = xlim[2], ytop = ylim[2],
        col = df1$rcol[i], border = df1$rcol[i]
      )

      # blank in part of main sector
      circos.rect(
        xleft = xlim[1], ybottom = ylim[1], xright = xlim[2] - rowSums(m)[i], ytop = ylim[1] + 0.3,
        col = "white", border = "white"
      )

      # white line all the way around
      circos.rect(xleft = xlim[1], ybottom = 0.3, xright = xlim[2], ytop = 0.32, col = "white", border = "white")

      # plot axis
      # circos.axis(labels.cex=0.6, direction = "outside", major.at=seq(from=0,to=floor(df1$xmax)[i],by=5),
      #            minor.ticks=1, labels.away.percentage = 0.15)
    }
  )
  ### Plot links (inner part)
  ### Add sum values to df1, marking the x-position of the first links
  ### out (sum1) and in (sum2). Updated for further links in loop below.
  df1$sum1 <- colSums(m)
  df1$sum2 <- numeric(n)

  ### Create a data.frame of the flow matrix sorted by flow size, to allow largest flow plotted first
  df2 <- cbind(as.data.frame(m), orig = rownames(m), stringsAsFactors = FALSE)
  df2 <- reshape(df2,
    idvar = "orig", varying = list(1:n), direction = "long",
    timevar = "dest", time = rownames(m), v.names = "m"
  )
  df2 <- arrange(df2, desc(m))

  ### Keep only the largest flows to avoid clutter
  df2 <- subset(df2, m > quantile(m, 0.6))
  
  ### Plot links
  for (k in 1:nrow(df2)) {
    # i,j reference of flow matrix
    i <- match(df2$orig[k], df1$ID)
    j <- match(df2$dest[k], df1$ID)

    # plot link
    circos.link(
      sector.index1 = df1$ID[i], point1 = c(df1$sum1[i], df1$sum1[i] + abs(m[i, j])),
      sector.index2 = df1$ID[j], point2 = c(df1$sum2[j], df1$sum2[j] + abs(m[i, j])),
      col = df1$lcol[i]
    )

    # update sum1 and sum2 for use when plotting the next link
    df1$sum1[i] <- df1$sum1[i] + abs(m[i, j])
    df1$sum2[j] <- df1$sum2[j] + abs(m[i, j])
  }
}
