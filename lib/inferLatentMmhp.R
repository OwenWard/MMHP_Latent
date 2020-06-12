myViterbi <- function(events, param, termination = NULL, if.print = FALSE){
  ## initialize
  if(is.null(termination)){
    termination <- tail(events, 1)+0.01
  }
  N <- length(events)
  interevent <- diff(c(0,events,termination))
  
  zt_v <- rep(NA,N)
  R <- rep(NA,N)
  
  m_matrix <- matrix(0,nrow=N,ncol=2)
  b_matrix <- matrix(0,nrow=N,ncol=2)
  m_termination <- rep(0,2)
  b_termination <- rep(0,2)
  
  probs_1_matrix <- matrix(0,nrow=N+1,ncol=2) ## Probability vector for transition to state 1
  probs_2_matrix <- matrix(0,nrow=N+1,ncol=2) ## Probability vector for transition to state 0
  int_1 <- matrix(0,nrow=N+1,ncol=2)
  int_2 <- matrix(0,nrow=N+1,ncol=2)
  logp1_1 <- rep(0,N+1)
  logp1_2 <- rep(0,N+1)
  logp2_1 <- rep(0,N+1)
  logp2_2 <- rep(0,N+1)
  
  ## --- Calculation for log probability of Markov transition logP_ij(t)
  for(n in 1:(N+1)){
    probs_1_matrix[n,1] <- log(param$q2/(param$q1+param$q2)+param$q1/(param$q1+param$q2)*exp(-(param$q1+param$q2)*interevent[n])) # 1->1
    probs_2_matrix[n,2] <- log(param$q1/(param$q1+param$q2)+param$q2/(param$q1+param$q2)*exp(-(param$q1+param$q2)*interevent[n])) # 2->2
    probs_1_matrix[n,2] <- log(1-exp(probs_2_matrix[n,2])) # 2->1
    probs_2_matrix[n,1] <- log(1-exp(probs_1_matrix[n,1])) # 1->2
  }
  ## --- R for Hawkes
  R[1] = 0;
  for(n in 2:(N+1)){
    R[n] = exp(-param$beta*interevent[n])*(R[n-1]+1)
  }
  ## --- Integration of lambda
  for(n in 1:(N+1)){
    K0 = exp(-(param$q1+param$q2)*interevent[n])
    K1 = (1-exp(-(param$q1+param$q2)*interevent[n]))/(param$q1+param$q2)
    K2 = (1-exp(-(param$q1+param$q2)*interevent[n]))/(param$q1+param$q2)
    K3 = R[n]*(exp(param$beta*interevent[n])-1)/param$beta
    K4 = R[n]*(1-exp(-(param$beta+param$q1+param$q2)*interevent[n]))*exp(param$beta*interevent[n])/(param$beta+param$q1+param$q2)
    K5 = R[n]*(1-exp(-(param$q1+param$q2-param$beta)*interevent[n]))/(param$q1+param$q2-param$beta)
    int_1[n,1] = ((param$q2^2*param$lambda1+param$q2*param$q1*param$lambda0)*interevent[n] +
                    (param$q1^2*param$lambda1+param$q2*param$q1*param$lambda0)*K0*interevent[n] +
                    (param$lambda1-param$lambda0)*param$q2*param$q1*K1 + (param$lambda1-param$lambda0)*param$q2*param$q1*K2 +
                    param$alpha*K3*(param$q2^2+param$q1^2*K0) +
                    param$alpha*param$q1*param$q2*K4 + param$alpha*param$q1*param$q2*K5)/(param$q1+param$q2)^2/exp(probs_1_matrix[n,1])  #1->1
    int_1[n,2] = ((param$q2^2*param$lambda1+param$lambda0*param$q1*param$q2)*interevent[n] -
                    (param$lambda1*param$q1*param$q2+param$lambda0*param$q2^2)*K0*interevent[n] +
                    (param$lambda0-param$lambda1)*param$q2^2*K1 + (param$lambda1-param$lambda0)*param$q1*param$q2*K2 +
                    param$alpha*param$q2*K3*(param$q2-param$q1*K0) -
                    param$alpha*param$q2^2*K4 + param$alpha*param$q1*param$q2*K5)/(param$q1+param$q2)^2/exp(probs_1_matrix[n,2]) #2->1
    int_2[n,1] = ((param$q1*param$q2*param$lambda1+param$q1^2*param$lambda0)*interevent[n] -
                    (param$q1^2*param$lambda1+param$q1*param$q2*param$lambda0)*K0*interevent[n] +
                    (param$lambda1-param$lambda0)*param$q1^2*K1 + param$q1*param$q2*(param$lambda0-param$lambda1)*K2 +
                    param$alpha*param$q1*K3*(param$q2-param$q1*K0) +
                    param$alpha*param$q1^2*K4 - param$alpha*param$q2*param$q1*K5)/(param$q1+param$q2)^2/exp(probs_2_matrix[n,1])  #1->2
    int_2[n,2] = ((param$q1*param$q2*param$lambda1+param$lambda0*param$q1^2)*interevent[n] +
                    (param$q1*param$q2*param$lambda1+param$lambda0*param$q2^2)*K0*interevent[n] +
                    (param$lambda0-param$lambda1)*param$q1*param$q2*K1 + (param$lambda0-param$lambda1)*param$q1*param$q2*K2 +
                    param$alpha*param$q1*param$q2*K3*(1+K0) -
                    param$alpha*param$q1*param$q2*K4 - param$alpha*param$q1*param$q2*K5)/(param$q1+param$q2)^2/exp(probs_2_matrix[n,2])  #2->2
  }
  
  ## --- consider n = 1
  #------------ transit to z(t_n) = 1
  logp1_1[1] <- probs_1_matrix[1,1] - int_1[1,1] # 1->1
  logp2_1[1] <- probs_1_matrix[1,2] - int_1[1,2] # 2->1
  if(logp1_1[1]>logp2_1[1]){
    m_matrix[1,1] <- logp1_1[1] + log(param$lambda1)
    b_matrix[1,1] <- 1
  }else{
    m_matrix[1,1] <- logp2_1[1] + log(param$lambda1)
    b_matrix[1,1] <- 2
  }
  
  #------------ transit to z(t_1) = 2
  logp1_2[1] <- probs_2_matrix[1,1] - int_2[1,1] # 1->2
  logp2_2[1] <- probs_2_matrix[1,2] - int_2[1,2] # 2->2
  if(logp1_2[1]>logp2_2[1]){
    m_matrix[1,2] <- logp1_2[1] + log(param$lambda0)
    b_matrix[1,2] <- 1
  }else{
    m_matrix[1,2] <- logp2_2[1] + log(param$lambda0)
    b_matrix[1,2] <- 2
  }
  
  if(N>1){
    ## interate n from 2 to N
    for(n in 2:N){
      #------------ transit to z(t_n) = 1
      logp1_1[n] <- m_matrix[n-1,1] + probs_1_matrix[n,1] - int_1[n,1] # 1->1
      logp2_1[n] <- m_matrix[n-1,2] + probs_1_matrix[n,2] - int_1[n,2] # 2->1
      if(logp1_1[n]>logp2_1[n]){
        m_matrix[n,1] <- logp1_1[n] + log(param$lambda1+param$alpha*R[n]) 
        b_matrix[n,1] <- 1
      }else{
        m_matrix[n,1] <- logp2_1[n] + log(param$lambda1+param$alpha*R[n]) 
        b_matrix[n,1] <- 2
      }
      
      #------------ transit to z(t_n) = 2
      logp1_2[n] <- m_matrix[n-1,1] + probs_2_matrix[n,1] - int_2[n,1] # 1->2
      logp2_2[n] <- m_matrix[n-1,2] + probs_2_matrix[n,2] - int_2[n,2] # 2->2
      if(logp1_2[n]>logp2_2[n]){
        m_matrix[n,2] <- logp1_2[n] + log(param$lambda0) 
        b_matrix[n,2] <- 1
      }else{
        m_matrix[n,2] <- logp2_2[n] + log(param$lambda0) 
        b_matrix[n,2] <- 2
      }
    }
  }
  
  ## termination
  #------------ transit to z(T) = 1
  logp1_1[N+1] <- m_matrix[N,1] + probs_1_matrix[N+1,1] - int_1[N+1,1] # 1->1
  logp2_1[N+1] <- m_matrix[N,2] + probs_1_matrix[N+1,2] - int_1[N+1,2] # 2->1
  if(logp1_1[N+1]>logp2_1[N+1]){
    m_termination[1] <- logp1_1[N+1]
    b_termination[1] <- 1
  }else{
    m_termination[1] <- logp2_1[N+1]
    b_termination[1] <- 2
  }
  
  #------------ transit to z(T) = 2
  logp1_2[N+1] <- m_matrix[N,1] + probs_2_matrix[N+1,1] - int_2[N+1,1] # 1->2
  logp2_2[N+1] <- m_matrix[N,2] + probs_2_matrix[N+1,2] - int_2[N+1,2] # 2->2
  if(logp1_2[N+1]>logp2_2[N+1]){
    m_termination[2] <- logp1_2[N+1]
    b_termination[2] <- 1
  }else{
    m_termination[2] <- logp2_2[N+1]
    b_termination[2] <- 2
  }
  
  # ----calulate zt_v for global decoding
  if(m_termination[1]>m_termination[2]){
    termination_state <- 1
  }else{
    termination_state <- 2
  }
  
  zt_v[N] <-  b_termination[termination_state]
  
  if(N > 1){
    for(k in 1:(N-1)){
      zt_v[N-k] = b_matrix[N-k+1,zt_v[N-k+1]];
    }
  }
  
  if(if.print){
    print_matrix <- rbind(events,t(m_matrix),t(b_matrix),zt_v)
    write.csv(print_matrix, "test.csv")
  }
  
  return(list(zt_v=zt_v, initial_state=b_matrix[1,zt_v[1]], termination_state=termination_state))
}

myViterbiWithInitial <- function(events, param, initial.p = 0.5, termination = NULL, if.print = FALSE){
  ## initialize
  if(is.null(termination)){
    termination <- tail(events, 1)
  }
  N <- length(events)
  interevent <- diff(c(0,events,termination))
  
  zt_v <- rep(NA,N)
  R <- rep(NA,N)
  
  m_matrix <- matrix(0,nrow=N,ncol=2)
  b_matrix <- matrix(0,nrow=N,ncol=2)
  m_termination <- rep(0,2)
  b_termination <- rep(0,2)
  
  probs_1_matrix <- matrix(0,nrow=N+1,ncol=2) ## Probability vector for transition to state 1
  probs_2_matrix <- matrix(0,nrow=N+1,ncol=2) ## Probability vector for transition to state 0
  int_1 <- matrix(0,nrow=N+1,ncol=2)
  int_2 <- matrix(0,nrow=N+1,ncol=2)
  logp1_1 <- rep(0,N+1)
  logp1_2 <- rep(0,N+1)
  logp2_1 <- rep(0,N+1)
  logp2_2 <- rep(0,N+1)
  
  ## --- Calculation for log probability of Markov transition logP_ij(t)
  for(n in 1:(N+1)){
    probs_1_matrix[n,1] <- log(param$q2/(param$q1+param$q2)+param$q1/(param$q1+param$q2)*exp(-(param$q1+param$q2)*interevent[n])) # 1->1
    probs_2_matrix[n,2] <- log(param$q1/(param$q1+param$q2)+param$q2/(param$q1+param$q2)*exp(-(param$q1+param$q2)*interevent[n])) # 2->2
    probs_1_matrix[n,2] <- log(1-exp(probs_2_matrix[n,2])) # 2->1
    probs_2_matrix[n,1] <- log(1-exp(probs_1_matrix[n,1])) # 1->2
  }
  ## --- R for Hawkes
  R[1] = 0;
  for(n in 2:(N+1)){
    R[n] = exp(-param$beta*interevent[n])*(R[n-1]+1)
  }
  ## --- Integration of lambda
  for(n in 1:(N+1)){
    K0 = exp(-(param$q1+param$q2)*interevent[n])
    K1 = (1-exp(-(param$q1+param$q2)*interevent[n]))/(param$q1+param$q2)
    K2 = (1-exp(-(param$q1+param$q2)*interevent[n]))/(param$q1+param$q2)
    K3 = R[n]*(exp(param$beta*interevent[n])-1)/param$beta
    K4 = R[n]*(1-exp(-(param$beta+param$q1+param$q2)*interevent[n]))*exp(param$beta*interevent[n])/(param$beta+param$q1+param$q2)
    K5 = R[n]*(1-exp(-(param$q1+param$q2-param$beta)*interevent[n]))/(param$q1+param$q2-param$beta)
    int_1[n,1] = ((param$q2^2*param$lambda1+param$q2*param$q1*param$lambda0)*interevent[n] +
                    (param$q1^2*param$lambda1+param$q2*param$q1*param$lambda0)*K0*interevent[n] +
                    (param$lambda1-param$lambda0)*param$q2*param$q1*K1 + (param$lambda1-param$lambda0)*param$q2*param$q1*K2 +
                    param$alpha*K3*(param$q2^2+param$q1^2*K0) +
                    param$alpha*param$q1*param$q2*K4 + param$alpha*param$q1*param$q2*K5)/(param$q1+param$q2)^2/exp(probs_1_matrix[n,1])  #1->1
    int_1[n,2] = ((param$q2^2*param$lambda1+param$lambda0*param$q1*param$q2)*interevent[n] -
                    (param$lambda1*param$q1*param$q2+param$lambda0*param$q2^2)*K0*interevent[n] +
                    (param$lambda0-param$lambda1)*param$q2^2*K1 + (param$lambda1-param$lambda0)*param$q1*param$q2*K2 +
                    param$alpha*param$q2*K3*(param$q2-param$q1*K0) -
                    param$alpha*param$q2^2*K4 + param$alpha*param$q1*param$q2*K5)/(param$q1+param$q2)^2/exp(probs_1_matrix[n,2]) #2->1
    int_2[n,1] = ((param$q1*param$q2*param$lambda1+param$q1^2*param$lambda0)*interevent[n] -
                    (param$q1^2*param$lambda1+param$q1*param$q2*param$lambda0)*K0*interevent[n] +
                    (param$lambda1-param$lambda0)*param$q1^2*K1 + param$q1*param$q2*(param$lambda0-param$lambda1)*K2 +
                    param$alpha*param$q1*K3*(param$q2-param$q1*K0) +
                    param$alpha*param$q1^2*K4 - param$alpha*param$q2*param$q1*K5)/(param$q1+param$q2)^2/exp(probs_2_matrix[n,1])  #1->2
    int_2[n,2] = ((param$q1*param$q2*param$lambda1+param$lambda0*param$q1^2)*interevent[n] +
                    (param$q1*param$q2*param$lambda1+param$lambda0*param$q2^2)*K0*interevent[n] +
                    (param$lambda0-param$lambda1)*param$q1*param$q2*K1 + (param$lambda0-param$lambda1)*param$q1*param$q2*K2 +
                    param$alpha*param$q1*param$q2*K3*(1+K0) -
                    param$alpha*param$q1*param$q2*K4 - param$alpha*param$q1*param$q2*K5)/(param$q1+param$q2)^2/exp(probs_2_matrix[n,2])  #2->2
  }
  
  ## --- consider n = 1
  #------------ transit to z(t_n) = 1
  logp1_1[1] <- log(initial.p) + probs_1_matrix[1,1] - int_1[1,1] # 1->1
  logp2_1[1] <- log(1-initial.p) + probs_1_matrix[1,2] - int_1[1,2] # 2->1
  if(logp1_1[1]>logp2_1[1]){
    m_matrix[1,1] <- logp1_1[1] + log(param$lambda1)
    b_matrix[1,1] <- 1
  }else{
    m_matrix[1,1] <- logp2_1[1] + log(param$lambda1)
    b_matrix[1,1] <- 2
  }
  
  #------------ transit to z(t_1) = 2
  logp1_2[1] <- log(initial.p) + probs_2_matrix[1,1] - int_2[1,1] # 1->2
  logp2_2[1] <- log(1-initial.p) + probs_2_matrix[1,2] - int_2[1,2] # 2->2
  if(logp1_2[1]>logp2_2[1]){
    m_matrix[1,2] <- logp1_2[1] + log(param$lambda0)
    b_matrix[1,2] <- 1
  }else{
    m_matrix[1,2] <- logp2_2[1] + log(param$lambda0)
    b_matrix[1,2] <- 2
  }
  
  if(N>1){
    ## interate n from 2 to N
    for(n in 2:N){
      #------------ transit to z(t_n) = 1
      logp1_1[n] <- m_matrix[n-1,1] + probs_1_matrix[n,1] - int_1[n,1] # 1->1
      logp2_1[n] <- m_matrix[n-1,2] + probs_1_matrix[n,2] - int_1[n,2] # 2->1
      if(logp1_1[n]>logp2_1[n]){
        m_matrix[n,1] <- logp1_1[n] + log(param$lambda1+param$alpha*R[n]) 
        b_matrix[n,1] <- 1
      }else{
        m_matrix[n,1] <- logp2_1[n] + log(param$lambda1+param$alpha*R[n]) 
        b_matrix[n,1] <- 2
      }
      
      #------------ transit to z(t_n) = 2
      logp1_2[n] <- m_matrix[n-1,1] + probs_2_matrix[n,1] - int_2[n,1] # 1->2
      logp2_2[n] <- m_matrix[n-1,2] + probs_2_matrix[n,2] - int_2[n,2] # 2->2
      if(logp1_2[n]>logp2_2[n]){
        m_matrix[n,2] <- logp1_2[n] + log(param$lambda0) 
        b_matrix[n,2] <- 1
      }else{
        m_matrix[n,2] <- logp2_2[n] + log(param$lambda0) 
        b_matrix[n,2] <- 2
      }
    }
  }
  
  ## termination
  #------------ transit to z(T) = 1
  logp1_1[N+1] <- m_matrix[N,1] + probs_1_matrix[N+1,1] - int_1[N+1,1] # 1->1
  logp2_1[N+1] <- m_matrix[N,2] + probs_1_matrix[N+1,2] - int_1[N+1,2] # 2->1
  if(logp1_1[N+1]>logp2_1[N+1]){
    m_termination[1] <- logp1_1[N+1]
    b_termination[1] <- 1
  }else{
    m_termination[1] <- logp2_1[N+1]
    b_termination[1] <- 2
  }
  
  #------------ transit to z(T) = 2
  logp1_2[N+1] <- m_matrix[N,1] + probs_2_matrix[N+1,1] - int_2[N+1,1] # 1->2
  logp2_2[N+1] <- m_matrix[N,2] + probs_2_matrix[N+1,2] - int_2[N+1,2] # 2->2
  if(logp1_2[N+1]>logp2_2[N+1]){
    m_termination[2] <- logp1_2[N+1]
    b_termination[2] <- 1
  }else{
    m_termination[2] <- logp2_2[N+1]
    b_termination[2] <- 2
  }
  
  # ----calulate zt_v for global decoding
  if(m_termination[1]>m_termination[2]){
    termination_state <- 1
  }else{
    termination_state <- 2
  }
  
  zt_v[N] <-  b_termination[termination_state]
  
  if(N > 1){
    for(k in 1:(N-1)){
      zt_v[N-k] = b_matrix[N-k+1,zt_v[N-k+1]];
    }
  }
  
  if(if.print){
    print_matrix <- rbind(events,t(m_matrix),t(b_matrix),zt_v)
    write.csv(print_matrix, "test.csv")
  }
  
  return(list(zt_v=zt_v, initial_state=b_matrix[1,zt_v[1]], termination_state=termination_state))
}

modifiedLatentTrajectory <- function(params=list(lambda0,lambda1,alpha,beta,q1,q2), interevent, zt, start=0){
  # input pars
  #       interevent: length N
  #       zt: inferred latent state, 1 or 2 , length N
  # output  z.hat: states of Markov process
  #         x.hat: time of each transition of Markov process
  lambda0 <- params$lambda0
  lambda1 <- params$lambda1
  alpha <- params$alpha
  beta <- params$beta
  q1 <- params$q1
  q2 <- params$q2
  N <- length(interevent)
  temp.t <- cumsum(interevent)
  if(length(unique(zt))==1){
    z.hat <- unique(zt)
    x.hat <- 0
  }else{
    z.hat <- rep(NA,sum(diff(zt)!=0)+1)
    x.hat <- rep(NA,sum(diff(zt)!=0)+1)
    z.hat[1] <- zt[1]
    x.hat[1] <- 0
    temp.count <- 1
    for(l in 2:N){
      if(zt[l]==1 & zt[l-1]==2){
        temp.count <- temp.count + 1
        #x.hat[temp.count] <- temp.t[l]
        #z.hat[temp.count] <- 1
        if((alpha<0)|(q2-q1<0)){
          x.hat[temp.count] <- temp.t[l-1]
          z.hat[temp.count] <- 1
        }else{
          temp.A <- alpha/beta*sum(exp(beta*temp.t[1:(l-1)]))
          temp.delta.t <- -temp.t[l-1]-log((q2-q1)/(temp.A*beta))/beta
          x.hat[temp.count] <- ifelse(temp.delta.t>0&(temp.delta.t<temp.t[l]-temp.t[l-1]),temp.t[l-1] + temp.delta.t,
                                      ifelse(temp.delta.t<=0,temp.t[l-1],temp.t[l]))
          z.hat[temp.count] <- 1
        }
      }
      if(zt[l]==2 & zt[l-1]==1){
        temp.count <- temp.count + 1
        x.hat[temp.count] <- temp.t[l-1]
        z.hat[temp.count] <- 2
      }
    }
  }
  return(list(x.hat=c(x.hat,tail(temp.t,1))+start,z.hat=c(z.hat,tail(z.hat,1))))
}

##Interpolate the latent state trajectory given estimated states of events

interpolateLatentTrajectory <- function(params=list(lambda0,lambda1,alpha,beta,q1,q2), 
                                        events, zt, initial.state=NULL, termination.time=NULL, termination.state=NULL,
                                        default.inactive = 2){
  # input params:
  #       events: length N (not include 0)
  #       zt: inferred latent state, 1 or 2 , length N
  #       initial.state: 1 or 2
  # output  z.hat: states of Markov process
  #         x.hat: time of each transition of Markov process
  
  interevent <- diff(c(0,events))
  N <- length(interevent)
  inactive_state <- setdiff(unique(zt),c(1))
  
  if(params$alpha < 0){
    params$alpha <- abs(params$alpha)
  }
  if(length(unique(zt))==1){
    ## no change at all
    z.hat <- rep(unique(zt),2)
    x.hat <- tail(events,1)
    
    ## termination state
    
    if(!is.null(termination.time)){
      if(is.null(termination.state)){
        if(unique(zt)==1){
          ## check whether state can change between last event and termination time, if change
          ## helper variables
          A.m <- cumsum(exp(params$beta*events)) #length = n; A=alpha/beta*A.m
          frequent.par <- params$q2-params$q1+params$lambda0-params$lambda1
          
          if(frequent.par<=0){
            x.hat[1] <- tail(events,1) # TODO
            z.hat[2] <- default.inactive
            x.hat[2] <- termination.time
            z.hat[3] <- default.inactive
          }else{
            if(is.finite(A.m[N])){
              l_0 <- params$alpha/params$beta*A.m[N]*exp(-params$beta*events[N])
              l_Delta <- frequent.par*(termination.time-events[N]) + params$alpha/params$beta*A.m[N]*exp(-params$beta*termination.time)
              if(l_0>l_Delta){
                x.hat[1] <- tail(events,1) 
                z.hat[2] <- default.inactive
                x.hat[2] <- termination.time
                z.hat[3] <- default.inactive
              }
            }
          }
        }
      }else{
        x.hat[1] <- tail(events,1) # TODO
        z.hat[2] <- termination.state
        x.hat[2] <- termination.time
        z.hat[3] <- termination.state
      }
    }
    
    ## initial state
    if(!is.null(initial.state)){
      if(initial.state!=unique(zt)){
        ## check whether state can change between 0 and first event time, if change
        # initial.delta <- events[1]/2 #TODO
        # x.hat <- c(initial.delta,x.hat)
        # z.hat <- c(initial.state,z.hat)
        
        ## if not change
      }
    }
  }else{
    z.hat <- rep(NA,sum(diff(zt)!=0)+1)
    x.hat <- rep(NA,sum(diff(zt)!=0)+1)
    
    ## helper variables
    A.m <- cumsum(exp(params$beta*events)) #length = n; A=alpha/beta*A.m
    frequent.par <- params$q2-params$q1+params$lambda0-params$lambda1
    
    z.hat[1] <- zt[1]
    #x.hat[1] <- 0
    temp.count <- 1
    for(l in 2:N){
      if(zt[l]==1 & zt[l-1]==inactive_state){ #inactive change to active
        if(frequent.par<=0){
          x.hat[temp.count] <- events[l-1]
        }else{
          temp.delta <- 1/params$beta*log(frequent.par/(A.m[l-1]*params$alpha))
          if(events[l-1]+temp.delta>=0){
            x.hat[temp.count] <- events[l-1]
          }else if(events[l]+temp.delta<=0){
            x.hat[temp.count] <- events[l]
          }else{
            x.hat[temp.count] <- -temp.delta
          }
        }
        z.hat[temp.count+1] <- 1
        temp.count <- temp.count + 1
      }
      
      if(zt[l]==inactive_state & zt[l-1]==1){ #active change to inactive
        if(frequent.par<=0){
          x.hat[temp.count] <- events[l-1]
        }else{
          l_0 <- params$alpha/params$beta*A.m[l-1]*exp(-params$beta*events[l-1])
          l_Delta <- frequent.par*(events[l]-events[l-1]) + params$alpha/params$beta*A.m[l-1]*exp(-params$beta*events[l])
          if(l_0>l_Delta){
            x.hat[temp.count] <- events[l-1]
          }else{
            x.hat[temp.count] <- events[l]
          }
        }
        
        z.hat[temp.count+1] <- inactive_state
        temp.count <- temp.count + 1
      }
    }
    
    ## termination
    if(is.null(termination.time)){
      x.hat[temp.count] <- tail(events,1)
      z.hat[temp.count+1] <- z.hat[temp.count]
    }else{
      if(is.null(termination.state)){
        if(z.hat[temp.count]==1){
          ## check whether state can change between last event and termination time, if change
          if(frequent.par<=0){
            x.hat[temp.count] <- tail(events,1)
            z.hat[temp.count+1] <- inactive_state
            x.hat[temp.count+1] <- termination.time
            z.hat[temp.count+2] <- inactive_state
          }else{
            l_0 <- params$alpha/params$beta*A.m[N]*exp(-params$beta*events[N])
            l_Delta <- frequent.par*(termination.time-events[N]) + params$alpha/params$beta*A.m[N]*exp(-params$beta*termination.time)
            if(is.finite(A.m[N])){
              if(l_0>l_Delta){
                x.hat[temp.count] <- tail(events,1)
                z.hat[temp.count+1] <- inactive_state
                x.hat[temp.count+1] <- termination.time
                z.hat[temp.count+2] <- inactive_state
              }else{
                x.hat[temp.count] <- termination.time
                z.hat[temp.count+1] <- z.hat[temp.count]
              }
            }else{
              x.hat[temp.count] <- termination.time
              z.hat[temp.count+1] <- z.hat[temp.count]
            }
          }
        }else{
          x.hat[temp.count] <- termination.time
          z.hat[temp.count+1] <- z.hat[temp.count]
        }
      }else{
        x.hat[temp.count] <- termination.time
        z.hat[temp.count+1] <- termination.state
      }
    }
    
    ## initial
    if(!is.null(initial.state)){
      if(initial.state!=zt[1]){
        ## check whether state can change between 0 and first event time, if change
        # initial.delta <- events[1]/2 #TODO
        # x.hat <- c(initial.delta, x.hat)
        # z.hat <- c(initial.state, z.hat)
        
        ## if not change
      }
    }
  }
  return(list(x.hat=x.hat,z.hat=z.hat))
}

