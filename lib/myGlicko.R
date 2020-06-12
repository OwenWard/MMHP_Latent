my_glicko <- function(x, status=NULL, init=c(2200,300), gamma=0, cval=15, history=FALSE, sort=FALSE, ...)
{ 
  gammas <- rep(gamma, length.out = nrow(x)) 
  names(x) <- c("Month","White","Black","Score")
  
  play <- sort(unique(c(x$White,x$Black)))
  np <- length(play)
  x$White <- match(x$White, play)
  x$Black <- match(x$Black, play)
  
  if(!is.null(status)) {
    npadd <- play[!(play %in% status$Player)]
    zv <- rep(0, length(npadd))
    npstatus <- data.frame(Player = npadd, Rating = rep(init[1],length(npadd)), 
                           Deviation = rep(init[2],length(npadd)), Games = zv, Win = zv, Draw = zv, 
                           Loss = zv, Lag = zv)
    if(!("Games" %in% names(status))) status <- cbind(status, Games = 0)
    if(!("Win" %in% names(status))) status <- cbind(status, Win = 0)
    if(!("Draw" %in% names(status))) status <- cbind(status, Draw = 0)
    if(!("Loss" %in% names(status))) status <- cbind(status, Loss = 0)
    if(!("Lag" %in% names(status))) status <- cbind(status, Lag = 0)
    status <- rbind(status[,c("Player","Rating","Deviation","Games","Win","Draw","Loss","Lag")], npstatus)
    rinit <- status[,2]
    dinit <- status[,3]
    ngames <- status[,4]
    nwin <- status[,5]
    ndraw <- status[,6]
    nloss <- status[,7]
    nlag <- status[,8]
    names(rinit) <- names(dinit) <- names(ngames) <- status$Player
  }
  else {
    rinit <- rep(init[1], length.out=np)
    dinit <- rep(init[2], length.out=np)
    ngames <- nwin <- ndraw <- nloss <- rep(0, length.out=np)
    nlag <- rep(0,np)
    names(rinit) <- names(dinit) <- names(ngames) <- names(nlag) <- play
  }
  
  nm <- length(unique(x$Month))
  curplay <- match(play, names(rinit))
  orats <- rinit[-curplay] 
  odevs <- dinit[-curplay]^2
  ongames <- ngames[-curplay]
  onwin <- nwin[-curplay]
  ondraw <- ndraw[-curplay]
  onloss <- nloss[-curplay]
  olag <- nlag[-curplay]
  olag[ongames != 0] <- olag[ongames != 0] + nm
  crats <- rinit[curplay] 
  cdevs <- dinit[curplay]^2
  ngames <- ngames[curplay] 
  nwin <- nwin[curplay]
  ndraw <- ndraw[curplay]
  nloss <- nloss[curplay]
  nlag <- nlag[curplay]
  
  qv <- log(10)/400; qv2 <- qv^2; qip3 <- 3*(qv/pi)^2 
  gammas <- split(gammas, x$Month)
  x <- split(x, x$Month)
  players <- lapply(x, function(y) unique(c(y$White, y$Black)))
  if(history) {
    histry <- array(NA, dim=c(np,nm,4), dimnames=list(play,1:nm,c("Rating","Deviation","Games","Lag")))
  }
  
  for(i in 1:nm) {
    traini <- x[[i]]
    gammai <- gammas[[i]]
    nr <- nrow(traini)
    playi <- players[[i]]
    
    cdevs[playi] <- pmin(cdevs[playi] + (nlag[playi]+1)*(cval^2), 122500)
    gdevs <- 1/sqrt(1 + qip3*cdevs) 
    ngamesi <- tabulate(c(traini$White,traini$Black), np)
    dscore <- .C("glicko_c",
                 as.integer(np), as.integer(nr), as.integer(traini$White-1), as.integer(traini$Black-1),
                 as.double(traini$Score), as.double(crats), as.double(gdevs), as.double(gammai),
                 dscore = double(2*np))$dscore
    dval <- dscore[(np+1):(2*np)]; dscore <- dscore[1:np]
    cdevs <- 1/(1/cdevs + dval)
    crats <- crats + cdevs * qv * dscore
    
    trainiplw <- c(traini$White[traini$Score==1],traini$Black[traini$Score==0])
    trainipld <- c(traini$White[traini$Score==0.5],traini$Black[traini$Score==0.5])
    trainipll <- c(traini$White[traini$Score==0],traini$Black[traini$Score==1])
    ngames <- ngames + ngamesi
    nwin <- nwin + tabulate(trainiplw, np)
    ndraw <- ndraw + tabulate(trainipld, np)
    nloss <- nloss + tabulate(trainipll, np)
    nlag[ngames!=0] <- nlag[ngames!=0] + 1
    nlag[playi] <- 0
    
    if(history) {
      histry[,i,1] <- crats
      histry[,i,2] <- sqrt(cdevs)
      histry[,i,3] <- ngames
      histry[,i,4] <- nlag
    }
  }
  if(!history) histry <- NULL
  player <- suppressWarnings(as.numeric(names(c(crats,orats))))
  if (any(is.na(player))) player <- names(c(crats,orats))
  dfout <- data.frame(Player=player, Rating=c(crats,orats), Deviation=sqrt(c(cdevs,odevs)), 
                      Games=c(ngames,ongames), Win=c(nwin,onwin), Draw=c(ndraw,ondraw), Loss=c(nloss,onloss), 
                      Lag=c(nlag,olag),
                      stringsAsFactors = FALSE)
  if(sort) dfout <- dfout[order(dfout$Rating,decreasing=TRUE),] else dfout <- dfout[order(dfout$Player),]
  row.names(dfout) <- 1:nrow(dfout)
  
  lout <- list(ratings=dfout, history=histry, gamma=gamma, cval=cval, type = "Glicko")
  class(lout) <- "rating"
  lout
}

myPlotGlicko <- function(df, cval=3, mycolors=c("black", "grey", "orange", "red"), 
                         ltypes=c(1,2,3,1,2,3,1,2,3,1,2,3), thetitle="",  linewd=1, ylim1=1000,ylim2=3250,
                         ndays=1){
  
  df <- as.data.frame(df)
  
  robj <- my_glicko(df, cval=cval, history=T)
  
  x<-as.data.frame(unlist(robj$history))  
  z<-as.factor(df[,1])  #this is the df the glicko was conducted on
  n<-nlevels(z)
  
  x.ratings<-x[,1:n]
  x.deviations<-x[,(1+n):(n+n)]
  
  #longform the data
  x.ratingsmelt<-reshape2::melt(x.ratings)
  
  ids<-rownames(x.ratings)       #to make id column
  x.ratingsmelt$ids<-rep(ids, n)  #making id column
  
  l.ids<-length(ids)
  x.ratingsmelt$event<-rep(1:n, each=l.ids) 
  
  
  #add ranks
  xrn<-as.data.frame(x.ratings[n])
  colnames(xrn)<-c("finalrating")
  
  x.ratingsmelt$rank<-rank(-xrn$finalrating, ties.method="random")
  
  #make ids1 a factor with levels defined by rank
  x.ratingsmelt1 <- data.frame(ids=unique(x.ratingsmelt$ids),rank=unique(x.ratingsmelt$rank))
  x.ratingsmelt1 <- x.ratingsmelt1[order(x.ratingsmelt1$rank),]
  x.ratingsmelt$ids1 <- factor(x.ratingsmelt$ids,levels=x.ratingsmelt1$ids)
  
  #define color palette, multiple options (see below)
  colourCount <-length(unique(x.ratingsmelt$ids))
  getPalette = colorRampPalette(mycolors)
  
  
  ### now plot using ids1 instead of ids.
  p1<-ggplot(x.ratingsmelt, aes(x = event, y = value, col=ids1, linetype=ids1)) +
    scale_colour_manual(values = getPalette(colourCount)) +
    scale_linetype_manual(values=ltypes) +
    ylab("Glicko Rating") +
    xlab("Event") +
    ggtitle(thetitle)+
    ylim(ylim1,ylim2)+
    geom_line(lwd = linewd) +
    theme(plot.title = element_text(hjust = 0, vjust = 1, size = rel(1.7)), 
          panel.background = element_blank(), 
          plot.background = element_blank(), 
          panel.grid.major.y = element_line(color = "gray75",linetype = 'dotted'), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank(), 
          strip.background  = element_blank(),
          strip.text = element_text(size=rel(1.1)),
          text = element_text(color="gray20", size=10),
          axis.text = element_text(size=rel(1.0)),
          axis.text.x = element_text(color="gray20", size=rel(1.0)),
          axis.text.y = element_text(color="gray20", size=rel(1.0)),
          axis.title.x = element_text(size=rel(1.0), vjust=0),
          axis.title.y = element_text(size=rel(1.0), vjust=1),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")
  
  
  return(p1)
}
