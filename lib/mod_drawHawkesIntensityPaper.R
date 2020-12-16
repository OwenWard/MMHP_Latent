mod_drawHawkesIntensityPaper<-function(lambda0,alpha,beta,events,color=1,line.width=1,
                                   add = TRUE,
                                   yupper=10,
                                   title_name="",title.cex=9, y.ratio=0, 
                                   min.y=NULL, min.x=0, max.x=NULL, 
                                   annotate.time=FALSE, event.cex=1,
                                   box.type="l", box.col="gray40", title.line=0.5){
  # input object: parameters for Hawkes process, include lambda0, alpha, beta 
  #       events: vector of event happening time
  horizon <- tail(events,1)
  events <- events[events>0]
  N <- length(events)
  n <- N
  
  if(add==FALSE){
    plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,yupper*(1+y.ratio)),
         type="n",xlab="",ylab="",xaxt="n",yaxt="n",axes=FALSE)
    box(lty ='solid', col="gray40", lwd=0.7, bty=box.type)
    title(title_name,line=title.line,cex.main=title.cex)
    #left axes
    #axis(2,at=seq(0,yupper,5),labels=FALSE,cex.lab=2,lwd=0,lwd.ticks=1) 
    #text(-termination/50, y=seq(0,yupper,5), labels=seq(0,yupper,5), cex=1.6, srt=0, xpd=TRUE)
    #bottom axes
    if(annotate.time){
      axis(1,at=seq(0,ceiling(termination),length.out=5),labels=FALSE,cex.lab=2,lwd=0,lwd.ticks=1) #bottom axes
      text(seq(0,ceiling(termination),length.out=5), y=-y.ratio*6, cex=1.7, srt=0, xpd=TRUE,
           labels=round(seq(0,ceiling(termination),length.out=5),0))
    }
    # events
    points(events[-1],rep(min.y,n-1),cex=event.cex,pch=16,col="blue")
    # points(state_time,rep(lambda0,m),cex=1.3,pch=4,col="red",lwd=2)
  }
  
  if(N==1){ # only one event condition
    segments(x0=0, x1=events[1], y0=lambda0, col=color,lwd=line.width)
    segments(x0=events[1], y0=lambda0, y1=lambda0+alpha, col=color,lwd=line.width)
    if(horizon > tail(events,1)){
      lambda.n <- function(s) lambda0+alpha*sum(exp(-beta*(s-events[1])))
      new.lambda.n<-Vectorize(lambda.n)
      curve(new.lambda.n, from=tail(events,1), to=horizon, add=TRUE, col=color,lwd=line.width)
    }
    if( !is.null(title) ) title(main=title)
  }else{
    segments(x0=0, x1=events[1], y0=lambda0, col=color,lwd=line.width)
    segments(x0=events[1], y0=lambda0, y1=lambda0+alpha, col=color,lwd=line.width)
    
    for(i in 1:(N-1)){
      lambda.n <- function(s) lambda0+alpha*sum(exp(-beta*(rep(s,i)-events[1:i])))
      new.lambda.n<-Vectorize(lambda.n)
      curve(new.lambda.n, from=events[i], to=events[i+1], add=TRUE, col=color,lwd=line.width)
      segments(x0=events[i+1],y0=lambda.n(events[i+1]),y1=lambda.n(events[i+1])+alpha, col=color,lwd=line.width)
    }
  }
}
