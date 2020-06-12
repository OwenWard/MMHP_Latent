formMatrix <- function(my.fun, my.vector){
  return.matrix <- matrix(0,nrow=length(my.vector),ncol=length(my.vector))
  for(i in c(1:length(my.vector))){
    for(j in c(1:length(my.vector))[-i]){
      return.matrix[i,j] <- my.fun(my.vector[i],my.vector[j])
    }
  }
  return(return.matrix)
}

simulateLatentHP <- function(lambda0_matrix,
                             alpha_matrix,
                             beta_matrix,
                             horizon){
  N <- nrow(lambda0_matrix)
  start <- numeric(0)
  end <- numeric(0)
  day_hour <- numeric(0)
  
  for(i in c(1:N)){
    for(j in c(1:N)[-i]){
      sim_uni_hawkes_object <- list(lambda0=lambda0_matrix[i,j],
                                    alpha=alpha_matrix[i,j],
                                    beta=beta_matrix[i,j])
      sim_uni_hawkes_events <- uniHawkesSimulation(object = sim_uni_hawkes_object, 
                                                   horizon = horizon)
      
      if(sim_uni_hawkes_events[1]!=0){
        start <- c(start,rep(i,length(sim_uni_hawkes_events)))
        end <- c(end,rep(j,length(sim_uni_hawkes_events)))
        day_hour <- c(day_hour,sim_uni_hawkes_events)
      }
    }
  }
  
  return(list(start=start,end=end,day_hour=day_hour))
}

simulateLatentMMHP <- function(lambda0_matrix,
                               lambda1_matrix,
                               alpha_matrix,
                               beta_matrix,
                               q1_matrix,
                               q2_matrix,
                               horizon,
                               if.prefer.active = FALSE){
  N <- nrow(lambda0_matrix)
  start <- numeric(0)
  end <- numeric(0)
  day_hour <- numeric(0)
  mmhp_matrix <- matrix(rep(list(), N*N),nrow = N, ncol =N) 
  for(i in c(1:N)){
    for(j in c(1:N)[-i]){
      #print(paste(i,"->",j))
      if(alpha_matrix[i,j]>0){
        sim_mmhp_object <- list(lambda0=lambda0_matrix[i,j],
                                lambda1=lambda1_matrix[i,j],
                                alpha=alpha_matrix[i,j],
                                beta=beta_matrix[i,j],
                                q1=q1_matrix[i,j],
                                q2=q2_matrix[i,j])
        if(sim_mmhp_object$alpha>sim_mmhp_object$beta){
          sim.cur.result <- NULL
          while(is.null(sim.cur.result)){
            sim.cur.result <- tryCatch(expr = {withTimeout(simulate.mmhp(object = sim_mmhp_object, 
                                                                         termination=horizon,
                                                                         if.prefer.active=if.prefer.active),
                                                           timeout = 5)},
                                       TimeoutException = function(ex) cat("Timeout. Skipping.\n"))
          }
          temp_mmhp <- sim.cur.result
        }else{
          temp_mmhp <- simulate.mmhp(object = sim_mmhp_object, 
                                     termination=horizon,
                                     if.prefer.active=if.prefer.active)
        }
        mmhp_matrix[i,j][[1]] <- temp_mmhp
        sim_uni_hawkes_events <- temp_mmhp$tau[temp_mmhp$tau>0]
        if(length(sim_uni_hawkes_events)>0){
          if(sim_uni_hawkes_events[1]!=0){
            start <- c(start,rep(i,length(sim_uni_hawkes_events)))
            end <- c(end,rep(j,length(sim_uni_hawkes_events)))
            day_hour <- c(day_hour,sim_uni_hawkes_events)
          }
        }
      }
    }
  }
  
  return(list(start=start, end=end, day_hour=day_hour, mmhp_matrix=mmhp_matrix))
}

cleanSimulationData <- function(raw_data, cut_off = 3, N=12){
  start <- raw_data$start[order(raw_data$day_hour)]
  end <- raw_data$end[order(raw_data$day_hour)]
  day_hour <- raw_data$day_hour[order(raw_data$day_hour)]
  M <- length(day_hour)
  
  # ------------------------------------
  # pair_count
  indicator_each_pair<-matrix(rep(list(),N*N),nrow=N, ncol=N)
  N_count<-matrix(rep(0,N^2),ncol=N,nrow=N)
  
  for(i in 1:N){
    for(j in c(1:N)[-i]){
      indicator_each_pair[i,j][[1]]<-which(start==i&end==j)
      N_count[i,j]<-length(unlist(indicator_each_pair[i,j]))
    }
  }
  
  # ------------------------------------
  # time_matrix & max_interevent
  I_fit <- matrix(rep(1:N,each=N),nrow=N,ncol=N,byrow=TRUE)[which(N_count>=cut_off)]
  J_fit <- matrix(rep(1:N,each=N),nrow=N,ncol=N)[which(N_count>=cut_off)]
  time_matrix <- matrix(0,ncol=max(N_count),nrow=sum(N_count>=cut_off))
  event_matrix <- matrix(0,ncol=max(N_count),nrow=sum(N_count>=cut_off))
  max_interevent <- rep(0,sum(N_count>=cut_off))
  
  for(i in 1:sum(N_count>=cut_off)){
    temp_t <- c(0,day_hour[unlist(indicator_each_pair[I_fit[i],J_fit[i]])])
    event_matrix[i,c(1:(length(temp_t)-1))] <- day_hour[unlist(indicator_each_pair[I_fit[i],J_fit[i]])]
    time_matrix[i,c(1:(length(temp_t)-1))] <- temp_t[-1]-temp_t[-length(temp_t)]
    max_interevent[i] <- max(temp_t[-1][-1]-temp_t[-1][-(length(temp_t)-1)])
  }
  return(list(N_count = N_count, time_matrix = time_matrix, event_matrix=event_matrix,
              max_interevent = max_interevent,
              I_fit = I_fit, J_fit = J_fit, day_hour = day_hour[1:M], 
              indicator_each_pair = indicator_each_pair, M = M,
              start = start[1:M], end = end[1:M]))
}


cleanSimulationDataForNCount <- function(raw_data, N = 12, n.start=0, n.end=NULL){
  if(is.null(n.end)){
    n.end <- max(raw_data$day_hour)
  }
  start <- raw_data$start[order(raw_data$day_hour)]
  end <- raw_data$end[order(raw_data$day_hour)]
  day_hour <- raw_data$day_hour[order(raw_data$day_hour)]
  
  start <- start[day_hour>n.start&day_hour<=n.end]
  end <- end[day_hour>n.start&day_hour<=n.end]
  day_hour <- day_hour[day_hour>n.start&day_hour<=n.end]
  M <- length(day_hour)
  
  # ------------------------------------
  # pair_count
  time_each_pair<-matrix(rep(list(),N*N),nrow=N, ncol=N)
  N_count<-matrix(rep(0,N^2),ncol=N,nrow=N)
  
  for(i in 1:N){
    for(j in c(1:N)[-i]){
      time_each_pair[i,j][[1]]<-day_hour[start==i&end==j]
      N_count[i,j]<-length(unlist(time_each_pair[i,j]))
    }
  }
  
  return(list(N_count = N_count, 
              time_each_pair = time_each_pair, 
              M = M))
}

glickoScoreSimulation <- function(raw_sim_data, history_df, history_horizon){

  simulate.idx <- (raw_sim_data$day_hour > history_horizon)
  df <- data.frame(Actor=raw_sim_data$start[simulate.idx],Recipient=raw_sim_data$end[simulate.idx],
                   Timestamp=raw_sim_data$day_hour[simulate.idx], score=rep(1,sum(simulate.idx)))
  df <- df[order(df$Timestamp),]
  df <- rbind(history_df[, c("Actor","Recipient","score")],df[, c("Actor","Recipient","score")])
  df$event <- 1:nrow(df)
  glick.df <- df[, c("event","Actor","Recipient","score"), with = FALSE]
  gl <- my_glicko(glick.df, history=TRUE, cval=2)

  return(gl)
}

glickoScoreSimulationWithWindowsReturnAll <- function(predict.sim, history.df){
  all.actor <- numeric(0)
  all.recipient <- numeric(0)
  for(w in c(1:length(predict.sim))){
    cur.order <- order(predict.sim[w][[1]]$day_hour)
    all.actor <- c(all.actor,predict.sim[w][[1]]$start[cur.order])
    all.recipient <- c(all.recipient,predict.sim[w][[1]]$end[cur.order])
  }
  
  df <- data.frame(Actor=all.actor,Recipient=all.recipient,score=rep(1,length(all.actor)))
  df <- rbind(history.df[, c("Actor","Recipient","score")],df[, c("Actor","Recipient","score")])
  df$event <- 1:nrow(df)
  return(my_glicko(df[, c("event","Actor","Recipient","score")], history=TRUE, cval=2))
}

glickoScoreSimulationWithWindows <- function(predict.sim, to.predice.obs, history.df){
  return.score.row.no <- rep(0,6)
  all.actor <- numeric(0)
  all.recipient <- numeric(0)
  for(cur.day in c(16:21)){
    observe.windows <- which(to.predice.obs$day==cur.day)
    for(w in observe.windows){
      return.score.row.no[cur.day-15] <- return.score.row.no[cur.day-15] + length(predict.sim[w][[1]]$start)
      cur.order <- order(predict.sim[w][[1]]$day_hour)
      all.actor <- c(all.actor,predict.sim[w][[1]]$start[cur.order])
      all.recipient <- c(all.recipient,predict.sim[w][[1]]$end[cur.order])
    }
  }
  
  df <- data.frame(Actor=all.actor,Recipient=all.recipient,score=rep(1,length(all.actor)))
  df <- rbind(history.df[, c("Actor","Recipient","score")],df[, c("Actor","Recipient","score")])
  df$event <- 1:nrow(df)
  gl.df <- my_glicko(df[, c("event","Actor","Recipient","score")], history=TRUE, cval=2)
  return(t(gl.df$history[,cumsum(return.score.row.no)+nrow(history.df),1]))
}

glickoScoreDSNL <- function(lambda.d.matrix, to.predice.obs,history.df){
  return.score.row.no <- rep(0,6)
  all.actor <- numeric(0)
  all.recipient <- numeric(0)
  
  for(cur.day in c(16:21)){
    observe.time <- sum(to.predice.obs[to.predice.obs$day==cur.day,"observe.time"])
    total.matrix <- round(lambda.d.matrix[cur.day-15,,]*observe.time)
    cur.actor <- rep(c(1:12),rowSums(total.matrix))
    cur.recipient <- rep(rep(c(1:12),12),as.vector(total.matrix))
    cur.order <- sample(c(1:length(cur.actor)))
    all.actor <- c(all.actor,cur.actor[cur.order])
    all.recipient <- c(all.recipient,cur.recipient[cur.order])
    return.score.row.no[cur.day-15] <- length(cur.actor)
  }
  
  df <- data.frame(Actor=all.actor,Recipient=all.recipient,score=rep(1,length(all.actor)))
  df <- rbind(history.df[, c("Actor","Recipient","score")],df[, c("Actor","Recipient","score")])
  df$event <- 1:nrow(df)
  gl.df <- my_glicko(df[, c("event","Actor","Recipient","score")], history=TRUE, cval=2)
  return(t(gl.df$history[,cumsum(return.score.row.no)+nrow(history.df),1]))
}
