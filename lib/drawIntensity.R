################################################################################################
#This R file is source code to generate intensity figure for the latex file
#There are two functions in this source code: drawUniMMHPIntensityPaper and drawHPIntensityPaper
#Function 'drawHPIntensityPaper' is used to assist main function 'drawUniMMHPIntensityPaper'
################################################################################################

drawUniMMHPIntensityPaper<-function(object,simulation,yupper=10,add=FALSE,color=1,line.width=1,
                                    title_name="",title.cex=9, y.ratio=0, min.y=NULL, min.x=0, max.x=NULL, 
                                    annotate.time=FALSE, event.cex=1, box.type="l", box.col="gray40", title.line=0.5){
  # input object: the parameter list used for generating mmhp
  #       simulation: simulation result from simulate.mmhp 
  
  t <- simulation$tau
  state <- simulation$z
  state_time <- simulation$x
  if(tail(state_time,1)<tail(t,1)){
    state_time <- c(state_time,tail(t,1))
    state <- c(state,3-tail(state,1))
  }
  lambda0 <- object$lambda0
  lambda1 <- object$lambda1
  alpha <- object$alpha
  beta <- object$beta
  termination <- max(tail(simulation$x,1),tail(simulation$tau,1))
  n <- length(t)
  m <- length(state)
  if(is.null(min.y)){
    min.y <- -lambda0
  }
  if(is.null(max.x)){
    max.x <- termination
  }
  #plot(1,4,xlim=c(0.27,23),ylim=c(-0.5,9), type="n",xlab="Time",ylab="",xaxt="n",yaxt="n",cex.lab=1.8,bty="n")
  
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
    points(t[-1],rep(min.y,n-1),cex=event.cex,pch=ifelse(simulation$zt[-1]==1,16,1),col="blue")
    points(state_time,rep(lambda0,m),cex=1.3,pch=4,col="red",lwd=2)
  }
  
  if((m==2&state[1]==1&state_time[2]==tail(t,1))|(length(unique(state))==1&state[1]==1)){
    drawHawkesIntensityPaper(lambda0,alpha,beta,t,color=color,line.width=line.width)
  }else{
    for(i in 1:(m-1)){
      if(state[i]==1){
        hawkes_time <- t[t>=state_time[i]&t<=state_time[i+1]]
        if(i==1) hawkes_time <- hawkes_time[-1]
        history <- t[t<state_time[i]]
        drawHPIntensityPaper(lambda1,i,alpha,beta,state_time[i],state_time[i+1],history[-1],hawkes_time,color=color,line.width=line.width)
      }else{
        segments(x0=state_time[i],x1=state_time[i+1], y0=lambda0, lty=1,col=color,lwd=line.width)
      }
      if(i>1){
        segments(x0=state_time[i],y0=lambda0, y1=lambda1, lty=2, col=color,lwd=line.width)
      }
    }
  }
  # if(add==FALSE){
  #   legend(9,32,c("State 1/0 events","State change point","True intensity","Estimated MMHP intensity",
  #                  "Estimated MMHPSD intensity","Estimated MMPP intensity"),
  #          col = c(NA,"red","black",rgb(0,190,0,maxColorValue = 255),
  #                  rgb(221,160,221,maxColorValue = 255),rgb(169,169,169,maxColorValue = 255)),
  #          y.intersp=0.88,x.intersp=0.15,
  #          pch = c(NA,4,NA,NA,NA,NA), pt.cex = c(NA,1.5,NA,NA,NA,NA),
  #          lty = c(NA,NA,1,1,1,1), lwd = c(NA,2.5,5,5,5,5),cex=1.6)
  #   points(9.43,30.09,pch=1,cex=1.8,col="blue")
  #   points(9.2,30.09,pch=16,cex=1.9,col="blue")
  # }
}

drawHPIntensityPaper<-function(lambda0,i,alpha,beta,start,end,history,hawkes_time,color=1,line.width=1){
  n <- length(hawkes_time)
  m <- length(history)
  
  if(n==0){
    if(i==1){
      segments(x0=start,x1=end,y0=lambda0,lwd=line.width)
    } else{
      lambda.n <- function(s) lambda0+alpha*sum(exp(-beta*(rep(s,m)-history)))
      new.lambda.n<-Vectorize(lambda.n)
      segments(x0=start,y0=lambda0,y1=lambda.n(end),lty=2,col=color,lwd=line.width)
      curve(new.lambda.n, from=start, to=end, add=TRUE,col=color,lwd=line.width)
    }
  }else{
    if(i==1){
      segments(x0=start,x1=hawkes_time[1],y0=lambda0,col=color,lwd=line.width)
      segments(x0=hawkes_time[1],y0=lambda0,y1=lambda0+alpha,col=color,lwd=line.width)
    } else{
      lambda.n <- function(s) lambda0+alpha*sum(exp(-beta*(rep(s,m)-history)))
      new.lambda.n<-Vectorize(lambda.n)
      segments(x0=start,y0=lambda0,y1=lambda.n(start),lty=2,col=color,lwd=line.width)
      curve(new.lambda.n, from=start, to=hawkes_time[1], add=TRUE,col=color,lwd=line.width)
      segments(x0=hawkes_time[1],y0=lambda.n(hawkes_time[1]),y1=lambda.n(hawkes_time[1])+alpha,col=color,lwd=line.width)
    }
    if(n>1){
      for(j in 1:(n-1)){
        lambda.n <- function(s) lambda0+alpha*sum(exp(-beta*(rep(s,m+j)-c(history,hawkes_time[1:j]))))
        new.lambda.n<-Vectorize(lambda.n)
        curve(new.lambda.n, from=hawkes_time[j], to=hawkes_time[j+1], add=TRUE,col=color,lwd=line.width)
        segments(x0=hawkes_time[j+1],y0=lambda.n(hawkes_time[j+1]),y1=lambda.n(hawkes_time[j+1])+alpha,col=color,lwd=line.width)
      }
    }
    lambda.n <- function(s) lambda0+alpha*sum(exp(-beta*(rep(s,m+n)-c(history,hawkes_time[1:n]))))
    new.lambda.n<-Vectorize(lambda.n)
    curve(new.lambda.n, from=hawkes_time[n], to=end, add=TRUE,col=color,lwd=line.width)
    segments(x0=end,y0=lambda.n(end),y1=lambda0,lty=2,col=color,lwd=line.width)
  }
}

drawHawkesIntensityPaper<-function(lambda0,alpha,beta,events,color=1,line.width=1){
  # input object: parameters for Hawkes process, include lambda0, alpha, beta 
  #       events: vector of event happening time
  horizon <- tail(events,1)
  events <- events[events>0]
  N <- length(events)
  
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

fixInterpolateState <- function(latent.inter, termination){
  x.hat <- c(0,latent.inter$x.hat[latent.inter$x.hat<termination])
  z.hat <- latent.inter$z.hat[-length(latent.inter$z.hat)]
  return(list(x.hat=x.hat, z.hat=z.hat))
}
fixStateTransition <- function(temp.mmhp){
  last.state <- temp.mmhp$z[1]
  last.time <- 0
  new.x <- c(last.time)
  new.z <- c(last.state)
  last.event.time <- 0
  for(x.idx in c(1:(length(temp.mmhp$x)))){
    if(x.idx==length(temp.mmhp$x)){
      time.range <- c(temp.mmhp$x[x.idx],tail(temp.mmhp$tau,1))
    }else{
      time.range <- c(temp.mmhp$x[x.idx],temp.mmhp$x[x.idx+1])
    }
    
    time.idx <- temp.mmhp$tau>time.range[1]&temp.mmhp$tau<time.range[2]
    if(sum(time.idx)>0){
      if(temp.mmhp$z[x.idx]!=last.state){
        last.time <- temp.mmhp$x[x.idx]
        last.state <- temp.mmhp$z[x.idx]
        new.x <- c(new.x,last.time)
        new.z <- c(new.z,last.state)
      }else if((temp.mmhp$tau[time.idx][1]-last.event.time>8)&last.state==1){
        new.x <- c(new.x,runif(1,last.event.time,last.event.time+3),
                   runif(1,temp.mmhp$tau[time.idx][1]-2,temp.mmhp$tau[time.idx][1]))
        new.z <- c(new.z,2,1)
      }
      last.event.time <- tail(temp.mmhp$tau[time.idx],1)
    }
  }
  # if(temp.mmhp$z[length(temp.mmhp$z)]==last.state){
  #   new.x <- c(new.x[1:(length(new.x)-1)],tail(temp.mmhp$x,1))
  # }else{
  #   new.x <- c(new.x,tail(temp.mmhp$x,1))
  #   new.z <- c(new.z,tail(temp.mmhp$z,1))
  # }
  
  return(list(x=new.x,z=new.z,tau=temp.mmhp$tau,zt=temp.mmhp$zt,lambda.max=temp.mmhp$lambda.max))
}
