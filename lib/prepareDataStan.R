############################################################################
## Function prepareDataStan
## Input: current_cohort | cohort number
## Output: N_til | #directed pair
##         no_observations| total number of observation windows
##         I_fit |vector corresponsing to the actors in time_matrix
##         J_fit |vector corresponsing to the recipients in time_matrix
##         max_Nm | maximum of number of events for each pair each window
##         Nm | number of events for each pair
##         interevent_time_matrix | include termination time difference in the last entry
##         event_matrix | event times in each observation window
##         delta_window | length of non-observation period #discard
##         finishing_time | for each observation window, the total observation time
############################################################################
prepareDataStan <- function(current_cohort){
  clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
  return_df <- cleanObservationPeriod(current_cohort,
                                      full_data[[cohort_names[current_cohort]]],
                                      clean_data)
  
  unique_pairs_df <- return_df %>% 
    group_by(initiator, recipient) %>%
    dplyr::summarize(count=n(), observe=list(observe.id))
  unique_observe_win <- unique(return_df[,c("observe.id",
                                            "observe.time",
                                            "observe.start",
                                            "observe.end")])
  
  time_array <- array(0, dim = c(nrow(unique_pairs_df),
                              max(return_df$observe.id),
                              max(unlist(lapply(return_df$event.times,length)))+1))
  event_array <- array(0, dim = c(nrow(unique_pairs_df),
                               max(return_df$observe.id),
                               max(unlist(lapply(return_df$event.times,length)))))
  count_matrix <- matrix(0, nrow = nrow(unique_pairs_df),
                         ncol = max(return_df$observe.id))
  delta_window <- c(0,
                    unique_observe_win$observe.start[-1] - 
                      unique_observe_win$observe.end[-nrow(unique_observe_win)])
  finishing_time <- matrix(0,
                           nrow=nrow(unique_pairs_df),
                           ncol=max(return_df$observe.id))
  max_interevent <- rep(0,nrow(unique_pairs_df))
  
  for(pair in 1:nrow(unique_pairs_df)){
    current_initiator <- as.numeric(unique_pairs_df[pair,"initiator"])
    current_recipient <- as.numeric(unique_pairs_df[pair,"recipient"])
    current_window_vec <- unique_pairs_df$observe[[pair]]
    for(win in 1:max(return_df$observe.id)){
      row_indicator <- return_df$initiator==current_initiator&return_df$recipient==current_recipient&return_df$observe.id==win
      if(sum(row_indicator)==0){ # no events during observation window
        count_matrix[pair,win] <- 0
        time_array[pair,win,1] <- unique_observe_win[unique_observe_win$observe.id==win,"observe.time"]
        max_interevent[pair] <- max(max_interevent[pair],time_array[pair,win,1])
      }else{
        time_vec <- return_df[row_indicator,"event.times"][[1]]
        observe_period <- unique_observe_win[unique_observe_win$observe.id==win,"observe.time"]
        count_matrix[pair,win] <- length(time_vec)
        time_array[pair,win,(1:(length(time_vec)+1))] <- diff(c(0,time_vec,observe_period))
        event_array[pair,win,(1:length(time_vec))] <- time_vec
        max_interevent[pair] <- max(max_interevent[pair],time_array[pair,win,(1:(length(time_vec)+1))])
      }
    }
  }
  return(list(N_til=nrow(unique_pairs_df),
              no_observations=max(return_df$observe.id),
              I_fit=unique_pairs_df$initiator,
              J_fit=unique_pairs_df$recipient,
              max_Nm=max(unlist(lapply(return_df$event.times,length))),
              Nm=count_matrix,
              interevent_time_matrix=time_array,
              event_matrix=event_array,
              #delta_window=delta_window,
              finishing_time=unique_observe_win$observe.time))
}

############################################################################
## Function prepareDataStanTrain
## Input: current_cohort | cohort number
##        train_day | 15 by default
## Output: N_til | #directed pair
##         no_observations| total number of observation windows
##         I_fit |vector corresponsing to the actors in time_matrix
##         J_fit |vector corresponsing to the recipients in time_matrix
##         max_Nm | maximum of number of events for each pair each window
##         Nm | number of events for each pair
##         interevent_time_matrix | include termination time difference in the last entry
##         event_matrix | event times in each observation window
##         delta_window | length of non-observation period # discard
##         finishing_time | for each observation window, the total observation time
############################################################################
prepareDataStanTrain <- function(current_cohort, train_day = 15, first_day = 1){
  clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
  return_df <- cleanObservationPeriod(current_cohort, 
                                      full_data[[cohort_names[current_cohort]]],
                                      clean_data)
  return_df <- return_df[return_df$day<=train_day, ]
  return_df <- return_df[return_df$day>= first_day, ]
  unique_pairs_df <- return_df %>%
    group_by(initiator, recipient) %>%
    dplyr::summarize(count = n(), observe = list(observe.id))
  unique_observe_win <- unique(return_df[,c("observe.id","observe.time",
                                            "observe.start","observe.end")])
  
  time_array <- array(0,
                      dim = c(nrow(unique_pairs_df),
                        length(unique(return_df$observe.id)),
                        max(unlist(lapply(return_df$event.times,length))) + 1 ))
  event_array <- array(0, dim = 
                         c(nrow(unique_pairs_df),
                           length(unique(return_df$observe.id)),
                           max(unlist(lapply(return_df$event.times,length)))))
  count_matrix <- matrix(0, 
                         nrow = nrow(unique_pairs_df),
                         ncol = length(unique(return_df$observe.id)))
  delta_window <- c(0, unique_observe_win$observe.start[-1] - 
                      unique_observe_win$observe.end[-nrow(unique_observe_win)])
  finishing_time <- matrix(0,
                           nrow = nrow(unique_pairs_df),
                           ncol = length(unique(return_df$observe.id)))
  max_interevent <- rep(0, nrow(unique_pairs_df))
  
  for(pair in 1:nrow(unique_pairs_df)){
    current_initiator <- as.numeric(unique_pairs_df[pair,"initiator"])
    current_recipient <- as.numeric(unique_pairs_df[pair,"recipient"])
    current_window_vec <- unique_pairs_df$observe[[pair]]
    for(win in seq_along(unique(return_df$observe.id))){
      win_id <- c(unique(return_df$observe.id))[win]
      row_indicator <- return_df$initiator==current_initiator&
        return_df$recipient==current_recipient&return_df$observe.id==win_id
      if(sum(row_indicator)==0){ # no events during observation window
        count_matrix[pair,win] <- 0
        time_array[pair,win,1] <- unique_observe_win[
          unique_observe_win$observe.id==win_id, "observe.time"]
        max_interevent[pair] <- max(max_interevent[pair],
                                    time_array[pair, win, 1])
      }else{
        time_vec <- return_df[row_indicator,"event.times"][[1]]
        observe_period <- unique_observe_win[
          unique_observe_win$observe.id==win_id, "observe.time"]
        count_matrix[pair,win] <- length(time_vec)
        time_array[pair,win,(1:(length(time_vec)+1))] <- diff(c(0,
                                                                time_vec,
                                                                observe_period))
        event_array[pair,win,(1:length(time_vec))] <- time_vec
        max_interevent[pair] <- max(max_interevent[pair],
                                    time_array[pair, win, 
                                               (1:(length(time_vec)+1))])
      }
    }
  }
  return(list(N_til = nrow(unique_pairs_df),
              no_observations = length(unique(return_df$observe.id)),
              I_fit = unique_pairs_df$initiator,
              J_fit = unique_pairs_df$recipient,
              max_Nm = max(unlist(lapply(return_df$event.times,length))),
              Nm = count_matrix,
              interevent_time_matrix = time_array,
              event_matrix = event_array,
              #delta_window=delta_window,
              finishing_time = unique_observe_win$observe.time))
}
############################################################################
## Function prepareDataStanIMMHP
## Input: current_cohort | cohort number
## Output: N_til | #directed pair
##         no_observations| total number of observation windows
##         max_Nm | maximum of number of events for each pair each window
##         Nm | number of events for each pair
##         time_matrix | include termination time difference in the last entry
##         max_interevent | for each pair, the maximum of interevent time
############################################################################
prepareDataStanIMMHP <- function(current_cohort){
  clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
  return_df <- cleanObservationPeriod(current_cohort,
                                      full_data[[cohort_names[current_cohort]]],
                                      clean_data)
  
  unique_pairs_df <- return_df %>% group_by(initiator, recipient) %>% dplyr::summarize(count=n(), observe=list(observe.id))
  unique_observe_win <- unique(return_df[,c("observe.id","observe.time")])
  
  time_array <- array(0,dim=c(nrow(unique_pairs_df),max(return_df$observe.id),max(unlist(lapply(return_df$event.times,length)))+1))
  count_matrix <- matrix(0,nrow=nrow(unique_pairs_df),ncol=max(return_df$observe.id))
  max_interevent <- rep(0,nrow(unique_pairs_df))
  
  for(pair in 1:nrow(unique_pairs_df)){
    current_initiator <- as.numeric(unique_pairs_df[pair,"initiator"])
    current_recipient <- as.numeric(unique_pairs_df[pair,"recipient"])
    current_window_vec <- unique_pairs_df$observe[[pair]]
    for(win in 1:max(return_df$observe.id)){
      row_indicator <- return_df$initiator==current_initiator&return_df$recipient==current_recipient&return_df$observe.id==win
      if(sum(row_indicator)==0){ # no events during observation window
        count_matrix[pair,win] <- 0
        time_array[pair,win,1] <- unique_observe_win[unique_observe_win$observe.id==win,"observe.time"]
        max_interevent[pair] <- max(max_interevent[pair],time_array[pair,win,1])
      }else{
        time_vec <- return_df[row_indicator,"event.times"][[1]]
        observe_period <- unique_observe_win[unique_observe_win$observe.id==win,"observe.time"]
        count_matrix[pair,win] <- length(time_vec)
        time_array[pair,win,(1:(length(time_vec)+1))] <- diff(c(0,time_vec,observe_period))
        max_interevent[pair] <- max(max_interevent[pair],time_array[pair,win,(1:(length(time_vec)+1))])
      }
    }
  }
  return(list(N_til=nrow(unique_pairs_df),
              no_observations=max(return_df$observe.id),
              max_Nm=max(unlist(lapply(return_df$event.times,length))),
              Nm=count_matrix,
              time_matrix=time_array,
              max_interevent=max_interevent))
}

############################################################################
## Function prepareDataStanTrainIMMHP
## Input: current_cohort | cohort number
##        train_day | 15 by default
## Output: N_til | #directed pair
##         no_observations| total number of observation windows
##         max_Nm | maximum of number of events for each pair each window
##         Nm | number of events for each pair
##         time_matrix | include termination time difference in the last entry
##         max_interevent | for each pair, the maximum of interevent time
############################################################################
prepareDataStanTrainIMMHP <- function(current_cohort, train_day = 15){
  clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
  return_df <- cleanObservationPeriod(current_cohort, 
                                      full_data[[cohort_names[current_cohort]]],
                                      clean_data)
  return_df <- return_df[return_df$day<=train_day,]
  
  unique_pairs_df <- return_df %>% group_by(initiator, recipient) %>% dplyr::summarize(count=n(), observe=list(observe.id))
  unique_observe_win <- unique(return_df[,c("observe.id","observe.time")])
  
  time_array <- array(0,dim=c(nrow(unique_pairs_df),max(return_df$observe.id),max(unlist(lapply(return_df$event.times,length)))+1))
  count_matrix <- matrix(0,nrow=nrow(unique_pairs_df),ncol=max(return_df$observe.id))
  max_interevent <- rep(0,nrow(unique_pairs_df))
  
  for(pair in 1:nrow(unique_pairs_df)){
    current_initiator <- as.numeric(unique_pairs_df[pair,"initiator"])
    current_recipient <- as.numeric(unique_pairs_df[pair,"recipient"])
    current_window_vec <- unique_pairs_df$observe[[pair]]
    for(win in 1:max(return_df$observe.id)){
      row_indicator <- return_df$initiator==current_initiator&return_df$recipient==current_recipient&return_df$observe.id==win
      if(sum(row_indicator)==0){ # no events during observation window
        count_matrix[pair,win] <- 0
        time_array[pair,win,1] <- unique_observe_win[unique_observe_win$observe.id==win,"observe.time"]
        max_interevent[pair] <- max(max_interevent[pair],time_array[pair,win,1])
      }else{
        time_vec <- return_df[row_indicator,"event.times"][[1]]
        observe_period <- unique_observe_win[unique_observe_win$observe.id==win,"observe.time"]
        count_matrix[pair,win] <- length(time_vec)
        time_array[pair,win,(1:(length(time_vec)+1))] <- diff(c(0,time_vec,observe_period))
        max_interevent[pair] <- max(max_interevent[pair],time_array[pair,win,(1:(length(time_vec)+1))])
      }
    }
  }
  return(list(N_til=nrow(unique_pairs_df),
              no_observations=max(return_df$observe.id),
              max_Nm=max(unlist(lapply(return_df$event.times,length))),
              Nm=count_matrix,
              time_matrix=time_array,
              max_interevent=max_interevent))
}
