############################################################################
## Function cleanData.R
## Input: raw_data | data.frame for one cohort, loaded from mice.Rdata
## Output: pair_count | 12*12 matrix, accumulated interaction for directed pair
##         time_matrix|#directed pair * max(#interactions among all pairs) matrix,
##                    |Each row in the matrix corresponds the accumulated interaction 
##                    |time (in which the unobserved period is removed)
##         max_interevent |vector with # directed pair entries
##                        |Each entry represents the maximum of the interevent time
##         I_fit |vector corresponsing to the actors in time_matrix
##         J_fit |vector corresponsing to the recipients in time_matrix
############################################################################

cleanData <- function(raw_data, cut_off = 3, N = 12){
  M_til <- dim(raw_data)[1]
  
  # ------------------------------------
  # Initialize M (# evetns), start, end
  start <- rep(0, M_til)
  end <- rep(0, M_til)
  day_hour <- rep(0, M_til)
  day <- rep(0, M_til)
  new_observation <- rep(0, M_til)
  time_stamp <- rep(0, M_til)
  observation_change <- numeric(0)
  
  start[1] <- as.numeric(as.character(raw_data[2,"Actor"]))
  end[1] <- as.numeric(as.character(raw_data[2,"Recipient"]))
  day_hour[1] <- raw_data[2,"newtime"]
  day[1] <- raw_data[2,"day"]
  new_observation[1] <- 1
  time_stamp[1] <- (12+raw_data[2,"hour"])+raw_data[2,"minute"]/60+raw_data[2,"secs"]/60^2
  M <- 1
  
  for(i in 3:M_til){
    if(raw_data[i,"Actor"] == "Start" | raw_data[i,"Actor"] == "End" |
       (raw_data[i-1,"newtime"] == raw_data[i,"newtime"]
        & raw_data[i-1,"Actor"] == raw_data[i,"Actor"]
        & raw_data[i-1,"Recipient"] == raw_data[i,"Recipient"])){
      next
    }else if(raw_data[i-2,"Actor"] == "End" & raw_data[i-1,"Actor"] == "Start"){
      M <- M + 1
      start[M] <- as.numeric(as.character(raw_data[i,"Actor"]))
      end[M] <- as.numeric(as.character(raw_data[i,"Recipient"]))
      day[M] <- raw_data[i,"day"]
      new_observation[M] <- new_observation[M-1] + 1
      time_stamp[M] <- (12+raw_data[i,"hour"])+raw_data[i,"minute"]/60+raw_data[i,"secs"]/60^2
      observation_change <- c(observation_change, day_hour[M-1] + raw_data[i-2,"newtime"] - raw_data[i-3,"newtime"])
      day_hour[M] <- day_hour[M-1] + raw_data[i,"newtime"] + raw_data[i-2,"newtime"] - raw_data[i-3,"newtime"]
    }else{
      M <- M + 1
      start[M] <- as.numeric(as.character(raw_data[i,"Actor"]))
      end[M] <- as.numeric(as.character(raw_data[i,"Recipient"]))
      day[M] <- raw_data[i,"day"]
      new_observation[M] <- new_observation[M-1]
      time_stamp[M] <- (12+raw_data[i,"hour"])+raw_data[i,"minute"]/60+raw_data[i,"secs"]/60^2
      day_hour[M] <- day_hour[M-1] + raw_data[i,"newtime"] - raw_data[i-1,"newtime"]
    }
  }
  day_hour <- day_hour/60
  
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
  max_interevent <- rep(0,sum(N_count>=cut_off))
  
  for(i in 1:sum(N_count>=cut_off)){
    temp_t <- c(0,day_hour[unlist(indicator_each_pair[I_fit[i],J_fit[i]])])
    time_matrix[i,c(1:(length(temp_t)-1))] <- temp_t[-1]-temp_t[-length(temp_t)]
    if(length(temp_t)>2){
      max_interevent[i] <- max(temp_t[-1][-1]-temp_t[-1][-(length(temp_t)-1)])
    }else{
      max_interevent[i] <- 0
    }
  }
  return(list(N_count = N_count,
              time_matrix = time_matrix,
              max_interevent = max_interevent,
              I_fit = I_fit,
              J_fit = J_fit,
              day_hour = day_hour[1:M], 
              indicator_each_pair = indicator_each_pair, M = M,
              start = start[1:M],
              end = end[1:M],
              day = day[1:M],
              new_observation = new_observation[1:M],
              time_stamp = time_stamp[1:M],
              observation_change = observation_change/60))
}

############################################################################
## Function cleanObservationPeriod.R
## Input: current_cohort | cohort number
##        clean_data | outpur from cleanData function
## Output: return_df | for each pair each observation period, record the event times
############################################################################

cleanObservationPeriod <- function(current_cohort, raw_df, clean_data, mice_number=12){
  # raw_df <- full_data[[cohort_names[current_cohort]]]
  # pass this as argument instead, so it can actually be modified
  return_df <- data.frame(day = numeric(),
                          observe.id = numeric(),
                          initiator = numeric(),
                          recipient = numeric(),
                          observe.start = numeric(),
                          observe.end = numeric(),
                          observe.time = numeric(),
                          observe.start.trans = numeric(),
                          observe.end.trans = numeric(),
                          no.events = numeric(),
                          event.times = numeric(),
                          event.times.trans = numeric())
  start_line <- which(raw_df$Actor=="Start"&raw_df$newtime==0)
  end_line <- which(raw_df$Actor=="End")
  cur_row <- 1
  for(current_obs in 1:length(start_line)){
    temp_df <- raw_df[(start_line[current_obs]+1):(end_line[current_obs]-1),]
    current_day <- unique(temp_df$day)
    temp_unique_pairs <- unique(temp_df[(temp_df$Actor%in%c(1:mice_number))&(temp_df$Recipient%in%c(1:mice_number)),
                                        c("Actor","Recipient")])
    
    return_df[cur_row:(cur_row+nrow(temp_unique_pairs)-1),"day"] <- unique(temp_df$day) #same
    return_df[cur_row:(cur_row+nrow(temp_unique_pairs)-1),"observe.id"] <- current_obs #same
    return_df[cur_row:(cur_row+nrow(temp_unique_pairs)-1),"observe.start"] <- as.numeric(strptime(as.character(raw_df[start_line[current_obs],"Timestamp"]), "%m/%d/%Y %H:%M:%OS")- 
                                                                                           strptime(as.character(raw_df[1,"Timestamp"]), "%m/%d/%Y %H:%M:%OS"),units = "hours")
    return_df[cur_row:(cur_row+nrow(temp_unique_pairs)-1),"observe.end"] <- as.numeric(strptime(as.character(raw_df[end_line[current_obs],"Timestamp"]), "%m/%d/%Y %H:%M:%OS")- 
                                                                                         strptime(as.character(raw_df[1,"Timestamp"]), "%m/%d/%Y %H:%M:%OS"),units = "hours")
    return_df[cur_row:(cur_row+nrow(temp_unique_pairs)-1),"observe.time"] <- raw_df[end_line[current_obs],"newtime"]/60 #tocheck
    return_df[cur_row:(cur_row+nrow(temp_unique_pairs)-1),"observe.start.trans"] <- c(0,clean_data$observation_change)[current_obs]
    return_df[cur_row:(cur_row+nrow(temp_unique_pairs)-1),"observe.end.trans"] <- c(0,clean_data$observation_change)[current_obs]+raw_df[end_line[current_obs],"newtime"]/60
    
    for(pair in 1:nrow(temp_unique_pairs)){
      return_df[cur_row,"initiator"] <- as.numeric(as.character(temp_unique_pairs[pair,"Actor"]))
      return_df[cur_row,"recipient"] <- as.numeric(as.character(temp_unique_pairs[pair,"Recipient"]))
      return_df[cur_row,"event.times"][[1]] <- list(unique(temp_df[(temp_df$Actor==as.numeric(as.character(temp_unique_pairs[pair,"Actor"])))&
                                                                     (temp_df$Recipient==as.numeric(as.character(temp_unique_pairs[pair,"Recipient"]))),"newtime"]/60))
      return_df[cur_row,"no.events"] <- length(return_df[cur_row,"event.times"][[1]])
      return_df[cur_row,"event.times.trans"][[1]] <- list(return_df[cur_row,"event.times"][[1]] + c(0,clean_data$observation_change)[current_obs])
      cur_row <- cur_row + 1
    }
  }
  
  return(return_df)
}

cleanDataForTraining <- function(all_raw_data, cut_off = 3, N = 12, train_period = 15){
  raw_data <- all_raw_data[all_raw_data[,"day"] <= train_period,]
  M_til <- dim(raw_data)[1]
  
  # ------------------------------------
  # Initialize M (# evetns), start, end
  start <- rep(0, M_til)
  end <- rep(0, M_til)
  day_hour <- rep(0, M_til)
  day <- rep(0, M_til)
  new_observation <- rep(0, M_til)
  time_stamp <- rep(0, M_til)
  
  start[1] <- as.numeric(as.character(raw_data[2,"Actor"]))
  end[1] <- as.numeric(as.character(raw_data[2,"Recipient"]))
  day_hour[1] <- raw_data[2,"newtime"]
  day[1] <- raw_data[2,"day"]
  new_observation[1] <- 1
  time_stamp[1] <- (12+raw_data[2,"hour"])+raw_data[2,"minute"]/60+raw_data[2,"secs"]/60^2
  M <- 1
  
  for(i in 3:M_til){
    if(raw_data[i,"Actor"] == "Start" | raw_data[i,"Actor"] == "End" |
       (raw_data[i-1,"newtime"] == raw_data[i,"newtime"]
        & raw_data[i-1,"Actor"] == raw_data[i,"Actor"]
        & raw_data[i-1,"Recipient"] == raw_data[i,"Recipient"])){
      next
    }else if(raw_data[i-2,"Actor"] == "End" & raw_data[i-1,"Actor"] == "Start"){
      M <- M + 1
      start[M] <- as.numeric(as.character(raw_data[i,"Actor"]))
      end[M] <- as.numeric(as.character(raw_data[i,"Recipient"]))
      day[M] <- raw_data[i,"day"]
      new_observation[M] <- new_observation[M-1] + 1
      time_stamp[M] <- (12+raw_data[i,"hour"])+raw_data[i,"minute"]/60+raw_data[i,"secs"]/60^2
      day_hour[M] <- day_hour[M-1] + raw_data[i,"newtime"] + raw_data[i-2,"newtime"] - raw_data[i-3,"newtime"]
    }else{
      M <- M + 1
      start[M] <- as.numeric(as.character(raw_data[i,"Actor"]))
      end[M] <- as.numeric(as.character(raw_data[i,"Recipient"]))
      day[M] <- raw_data[i,"day"]
      new_observation[M] <- new_observation[M-1]
      time_stamp[M] <- (12+raw_data[i,"hour"])+raw_data[i,"minute"]/60+raw_data[i,"secs"]/60^2
      day_hour[M] <- day_hour[M-1] + raw_data[i,"newtime"] - raw_data[i-1,"newtime"]
    }
  }
  day_hour <- day_hour/60
  
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
  max_interevent <- rep(0,sum(N_count>=cut_off))
  
  for(i in 1:sum(N_count>=cut_off)){
    if(N_count[I_fit[i],J_fit[i]]==0){
      max_interevent[i] <- tail(day_hour[1:M],1)
    }else if(N_count[I_fit[i],J_fit[i]]==1){
      time_matrix[i,1] <- day_hour[unlist(indicator_each_pair[I_fit[i],J_fit[i]])]
      max_interevent[i] <- max(tail(day_hour[1:M],1)-day_hour[unlist(indicator_each_pair[I_fit[i],J_fit[i]])], day_hour[unlist(indicator_each_pair[I_fit[i],J_fit[i]])])
    }else{
      temp_t <- c(0,day_hour[unlist(indicator_each_pair[I_fit[i],J_fit[i]])])
      time_matrix[i,c(1:(length(temp_t)-1))] <- temp_t[-1]-temp_t[-length(temp_t)]
      max_interevent[i] <- max(temp_t[-1][-1]-temp_t[-1][-(length(temp_t)-1)])
    }
  }
  return(list(N_count = N_count, time_matrix = time_matrix, max_interevent = max_interevent,
              I_fit = I_fit, J_fit = J_fit, day_hour = day_hour[1:M], 
              indicator_each_pair = indicator_each_pair, M = M,
              start = start[1:M], end = end[1:M],
              day = day[1:M], new_observation = new_observation[1:M], time_stamp = time_stamp[1:M]))
}
