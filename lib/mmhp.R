################################################################################################
#This R file is source code to including functions related to Markov Modulated Hawkes Process (MMHP)
#1. simulated.mmhp: to simulate MMHP using thinning theorem
#2. simulate.hp: to simulate Hawkes process during state 1 (used in function 'simulated.mmhp')
#3. simulate.mmhp.given.state: to simulate MMHP given hidden state trajectory
#4. drawUniMMHPIntensity: plot intensity function for MMHP
#5. drawHPIntensity: a helper function for 'drawUniMMHPIntensity', to draw Hawkes process intensity
#6. mmhpNegLogLik: compute negative log likelihood for mmhp
#7. mmhpCompensator: compute vector of compensators for mmhp
#8. mmhpIntensityNumeric: Compute numerical value of intensity function
#9. mmhpIntensityAtEvents:Compute intensity evaluates at event times (output is an vector)
#10. mmhpViterbi: infer the probability of zt=1 for each event using Viterbi algorithm
#11. mmhpLocalState: infer the probability of zt=1 for each event using local decoding
#12. modifiedLatentTrajectory: interpolate the latent state trajectory given estimated states of events
################################################################################################

##Simulate Markov Modulated Hawkes Process (including all the history)
simulate.mmhp <- function(object, termination, nsim=1000, seed=NULL, if.prefer.active = FALSE, ...){
  # VALID FOR 2 STATES: 1: active, 2: not active
  # input q1,q2: transition probability
  #       delta: initial state probability; 
  #       lambda1, alpha, beta : parameters for Hawkes process
  #       lambda0: parameters for homogeneous Poisson process
  #       nsim: preset number of events
  # output  z: states of Markov process
  #         x: time of each transition of Markov process
  #         zt: state at each event
  #         tau: times of Poisson events
  if (!is.null(seed)) set.seed(seed)
  m <- 2
  #------------------------
  if(is.null(object$delta)){
    if(if.prefer.active){
      initial <- 1
    }else{
      initial <- sample(m, 1)
    }
  }else if (sum(object$delta)!=1){
    stop("Invalid delta")
  }
  
  if (any(object$delta==1)){
    initial <- (1:m)[as.logical(object$delta)]
  }else{
    initial <- sample(m, 1, prob=object$delta)
  }
  
  #------------------------
  Q <- matrix(c(-object$q1,object$q1,object$q2,-object$q2),ncol=m,byrow=TRUE)
  lambda0 <- object$lambda0
  lambda1 <- object$lambda1
  alpha <- object$alpha
  beta <- object$beta
  
  Pi <- diag(m) - diag(1/diag(Q)) %*% Q
  zt <- rep(NA, nsim+1)
  tau <- rep(NA, nsim+1)
  #------------------------ initialization for Markov process
  #    the length of x and z may be too short
  #    gets extended later if required
  x <- rep(NA, nsim*10)
  z <- rep(NA, nsim*10)
  z[1] <- zt[1] <- initial
  x[1] <- tau[1] <- 0
  last.t <- 0
  lambda.max <- 0
  i <- 1 #index for state
  j <- 2 #index for event
  #------------------------ initialization for Hawkes process
  
  while (j < nsim+2 & last.t < termination){
    i <- i+1
    #   extend x and z if too short
    if (i > length(x)){
      x <- c(x, rep(NA, nsim*10))
      z <- c(z, rep(NA, nsim*10))
    }
    #   sim time spent in Markov state y[i-1]
    z[i] <- sample(x=1:m, size=1, prob=Pi[(z[i-1]),])
    x[i] <- x[i-1] + rexp(1, rate=-Q[z[i-1], z[i-1]])
    t0 <- x[i-1]
    
    if(z[i-1]==1){
      #   sim times of Hawkes Poisson events
      if(x[i]>termination){
        x[i] <- termination
      }
      simulate.result <- simulate.hp(lambda1,alpha,beta,x[i-1],x[i],tau[1:(j-1)])
      hp <- simulate.result$t
      lambda.max <- ifelse(lambda.max>simulate.result$lambda.max,lambda.max,simulate.result$lambda.max)
      if(!hp[1]==0){
        tau[j:(j+length(hp)-1)] <- hp
        zt[j:(j+length(hp)-1)] <- z[i-1]
        j <- j + length(hp)
        last.t <- tail(hp,1)
      }
    }
    
    if(z[i-1]==2){
      while(j < nsim+2 & last.t < termination){
        #   sim times of Poisson events
        ti <- t0 + rexp(1, rate=lambda0)
        if (ti < x[i]){
          tau[j] <- t0 <- ti
          zt[j] <- z[i-1]
          j <- j + 1
          last.t <- ti
        }
        else break
      }
    }
  }
  
  final_i <- min(i,sum(x[1:i]<termination))
  final_nsim <- min(nsim,sum(tau<=termination,na.rm=TRUE))
  return(list(x=x[1:final_i],z=z[1:final_i],
              tau=tau[1:final_nsim],zt=zt[1:final_nsim],
              lambda.max=lambda.max))
}

##Simulate Hawkes process during active state (including all the histoty)
simulate.hp <- function(lambda0,alpha,beta,start,horizon,history){
  j0<-length(history)+1
  lambda.star<-ifelse(j0==2,lambda0,lambda0+alpha*sum(exp(-beta*(rep(start,j0-2)-history[2:(j0-1)]))))
  lambda.max <- lambda.star
  t<-numeric(10)
  n<-1
  U<-runif(1)
  s<--log(U)/lambda.star
  ti<-start+s
  repeat {
    if (ti > horizon){
      break
    } 
    
    lambda.star<-lambda.star+alpha
    t[n]<-ti
    if (length(t) < n+1) t <- c(t, numeric(10)) 
    
    repeat{
      U<-runif(1)
      s<-s-log(U)/lambda.star
      ti<-start+s
      lambda.s<-lambda0+alpha*sum(exp(-beta*c(rep(ti,n)-t[1:n],rep(ti,j0-1)-history[1:j0-1])))
      D<-runif(1)
      if(D<=lambda.s/lambda.star){
        lambda.star<-lambda.s
        lambda.max <- ifelse(lambda.max>lambda.star,lambda.max,lambda.star)
        break
      }
      lambda.star<-lambda.s
      lambda.max <- ifelse(lambda.max>lambda.star,lambda.max,lambda.star)
    }
    
    n <- n + 1
  }
  
  return(list(t=t[1:(n-1)],lambda.max=lambda.max))
}

hawkesIntensityNumeric <- function(params=list(lambda1,alpha,beta), t, time.vec){
  lambda1.t <- rep(0,length(time.vec))
  event.idx <- 1
  
  r <- 0
  for(i in c(1:length(time.vec))){
    current.t <- time.vec[i]
    if(event.idx < length(t)){
      if(current.t>t[event.idx+1]){
        event.idx <- event.idx + 1
        r <- exp(-params$beta*(t[event.idx]-t[event.idx-1]))*(1+r)
      }
    }
    
    if(current.t<=t[1]){
      lambda1.t[i]<-params$lambda1
    }else{
      lambda1.t[i]<-params$lambda1+params$alpha*exp(-params$beta*(current.t-t[event.idx]))*(1+r)
    }
  }
  
  return(lambda1.t)
}

mmhpIntensityNumeric <-function(params=list(lambda0,lambda1,alpha,beta,q1,q2), t, time.vec, latent.vec){
  ## latent.vec is vector with same length as time.vec, each entry is the probability at state 1
  lambda1.t <- hawkesIntensityNumeric(params=list(lambda1=params$lambda1,
                                                  alpha=params$alpha,
                                                  beta=params$beta), t, time.vec)
  lambda.t <- lambda1.t*latent.vec + params$lambda0*(1-latent.vec)
  return(lambda.t)
}


mmhpIntensityNumeric_win <-function(params=list(lambda0,lambda1,alpha,beta,q1,q2), t, time.vec, latent.vec){
  ## latent.vec is vector with same length as time.vec, each entry is the probability at state 1
  lambda1.t <- hawkesIntensityNumeric_win(params=list(lambda1=params$lambda1,
                                                      alpha=params$alpha,
                                                      beta=params$beta), t, time.vec)
  lambda.t <- lambda1.t*latent.vec + params$lambda0*(1-latent.vec)
  return(lambda.t)
}



##Compute numerical value of intensity function
mmhpTrueIntensityNumeric <- function(params=list(lambda0,lambda1,alpha,beta,q1,q2), t, 
                                     latent, time.vec=NULL){
  if(is.null(time.vec)){
    delta.t <- tail(t,1)/20000
    time.vec <- seq(delta.t,tail(t,1),delta.t)
  }
  lambda.t <- rep(params$lambda0,length(time.vec))
  x.idx <- 2
  event.idx <- 2
  if(tail(latent$x,1)<tail(t,1)){
    latent <- list(x=c(latent$x,tail(t,1)),
                   z=c(latent$z,3-tail(latent$z,1)))
  }
  r <- 0
  for(i in c(1:length(time.vec))){
    current.t <- time.vec[i]
    if(current.t>t[event.idx+1]){
      event.idx <- event.idx + 1
      r <- exp(-params$beta*(t[event.idx]-t[event.idx-1]))*(1+r)
    }
    if(current.t>latent$x[x.idx]){
      x.idx <- x.idx + 1
    }
    if(latent$z[x.idx-1]==1){
      if(current.t<=t[2]){
        lambda.t[i]<-params$lambda1
      }else{
        lambda.t[i]<-params$lambda1+params$alpha*exp(-params$beta*(current.t-t[event.idx]))*(1+r)
      }
    }
  }
  if(is.null(time.vec)){
    return(list(lambda.t=lambda.t,delta.t=delta.t,time.vec=time.vec))
  }else{
    return(list(lambda.t=lambda.t,time.vec=time.vec))
  }
}

##Compute intensity evaluates at event times (output is an vector)
mmhpIntensityAtEvents<-function(params=list(lambda0,lambda1,alpha,beta,q1,q2), t, latent_z){
  if(t[1]==0){
    t <- t[-1]
  }
  if(length(latent_z)==(length(t)+1)){
    latent_z <- latent_z[-1]
  }
  lambda.t <- rep(params$lambda0,length(t))
  r <- 0
  for(i in c(1:length(t))){
    if(i>1){
      r <- exp(-params$beta*(t[i]-t[i-1]))*(1+r)
    }
    if(latent_z[i]==1){
      lambda.t[i]<-params$lambda1+params$alpha*r
    }
  }
  return(lambda.t)
}

##Compute rescaled interevent time (input is from mmhpIntensityNumeric)
mmhpRescaledInterevent <- function(event.time, lambda.numeric.vec, time.segment.vec){
  ## event time not include 0
  Lambda_vec <- rep(0,length(event.time)-1)
  for(m in 1:(length(event.time)-1)){
    interevent_idx <- (time.segment.vec>=event.time[m]) & (time.segment.vec<event.time[m+1])
    Lambda_vec[m] <- sum(lambda.numeric.vec[interevent_idx])*(time.segment.vec[2]-time.segment.vec[1])
  }
  return(Lambda_vec)
}

##Compute vector of compensators for mmhp
mmhpCompensator<-function(params=list(lambda0,lambda1,alpha,beta,q1,q2), t, pzt, termination=NULL,if.pzt=TRUE){
  lambda0 <- params$lambda0
  lambda1 <- params$lambda1
  alpha <- params$alpha
  beta <- params$beta
  q1 <- params$q1
  q2 <- params$q2
  n <- length(t) - 1
  interevent <- t[-1]-t[-(n+1)]
  if(if.pzt){
    pzt <- mmhpLocalState(params=params, interevent)$pzt
  }
  
  ##compute compensator for Hawkes process
  Lambda <- rep(0,n)
  A<-0
  Lambda[1]<-lambda0*(interevent[1])
  for(i in 2:n){
    A<-1+exp(-beta*(interevent[i-1]))*A
    Lambda[i]<-lambda1*(interevent[i])+alpha/beta*(1-exp(-beta*interevent[i]))*A
  }
  
  Lambda_mixed <- Lambda*pzt + lambda0*interevent*(1-pzt)
  
  if(!is.null(termination)){
    A<-1+exp(-beta*(termination-tail(t,1)))*A
    last.lambda<-(lambda0*(termination-tail(t,1))+lambda1*(termination-tail(t,1))+alpha/beta*(1-exp(-beta*(termination-tail(t,1))))*A)/2
    return(c(Lambda_mixed,last.lambda))
  }else{
    return(Lambda_mixed)
  }
}

mmhpIntensityAtTime<-function(params=list(lambda0,lambda1,alpha,beta,q1,q2), events, current_time, latent_z, latent_x){
  if(is.null(latent_z)){
    return(params$lambda0)
  }
  if(current_time<latent_x[1]){
    current_state <- latent_z[1]
  }else{
    if(current_time>tail(latent_x,1)){
      current_state <- tail(latent_z,1)
    }else{
      current_state <- latent_z[-length(latent_x)][current_time>latent_x[-length(latent_x)]&current_time<=latent_x[-1]]
    }
  } 
  
  if(current_state==2){
    lambda.t <- params$lambda0
  }else{
    lambda.t <- uniHawkesIntensity(object=list(lambda0=params$lambda1,alpha=params$alpha,beta=params$beta), 
                                   events=events, 
                                   current_time=current_time)
  }
  return(lambda.t)
}