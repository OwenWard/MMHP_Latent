prepareData_permute_Stan <- function(current_cohort) {
  ### need to permute full_data[[cohort_names[current_cohort]]]
  ### and then go from there
  
  permuted <- permute_raw(current_cohort)
  clean_data <- cleanData(permuted)
  return_df <- cleanObservationPeriod(current_cohort, clean_data)
  
  unique_pairs_df <- return_df %>% group_by(initiator, recipient) %>% dplyr::summarize(count=n(), observe=list(observe.id))
  unique_observe_win <- unique(return_df[,c("observe.id","observe.time","observe.start","observe.end")])
  
  time_array <- array(0,dim=c(nrow(unique_pairs_df),max(return_df$observe.id),max(unlist(lapply(return_df$event.times,length)))+1))
  event_array <- array(0,dim=c(nrow(unique_pairs_df),max(return_df$observe.id),max(unlist(lapply(return_df$event.times,length)))))
  count_matrix <- matrix(0,nrow=nrow(unique_pairs_df),ncol=max(return_df$observe.id))
  delta_window <- c(0,unique_observe_win$observe.start[-1]-unique_observe_win$observe.end[-nrow(unique_observe_win)])
  finishing_time <- matrix(0,nrow=nrow(unique_pairs_df),ncol=max(return_df$observe.id))
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
              finishing_time=unique_observe_win$observe.time,
              permuted_data = permuted))
}


permute_raw <- function(current_cohort) {
  orig_raw <- full_data[[cohort_names[current_cohort]]]
  # then just permute the actor and recipient of these,
  # excluding the start and end rows
  start_rows <- which(orig_raw$Actor == "Start")  
  end_rows <- which(orig_raw$Actor == "End")
  
  exclude <- c(start_rows,end_rows)
  
  rows <- 1:nrow(orig_raw)
  shuffle_rows <- rows[-exclude]
  
  mice <- 1:12
  rand_actor <- sample(mice, size = length(shuffle_rows),replace = T)
  rand_recip <- rep(NA,length(shuffle_rows))
  for(i in 1:length(rand_actor)) {
    rand_recip[i] <- sample(mice[-rand_actor[i]],size = 1)
  }
  
  orig_raw$Actor[shuffle_rows] <- rand_actor
  orig_raw$Recipient[shuffle_rows] <- rand_recip
  
  return(orig_raw)
}



