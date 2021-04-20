# this file refits the complete I-MMHP model
# including prediction and Pearson residuals for each cohort


### code ###
## run this section if running on the cluster
source("/rigel/stats/users/ogw2103/code/MMHP/MMHP_Latent/run_scripts/cluster_setup.R")
source("run_scripts/cluster_setup.R")
### set cohort_id based on job num
jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
cohort_id <- jobid
#cohort_id <- 1
#####

save_data_path <- "output/common_rate/"

no_segments <- 500

library(rstan)
options(mc.cores = parallel::detectCores())


library(compete)
#library(RColorBrewer)
#library(Hmisc)
#library(wCorr)
#library(tidyverse)
library(dplyr)
#library(R.utils)
#library(fields)



source('lib/naiveRankHierarchy.R')
source('lib/expertRankHierarchy.R')
source('lib/cleanData.R') 
source('lib/prepareDataStan.R') 
source('lib/inferLatentMmhp.R')
source('lib/plotUtil.R') 
source('lib/mmhp.R')
source('lib/uniHawkes.R')
source('lib/simulatePrediction.R')
source('lib/myGlicko.R')
source("https://gist.githubusercontent.com/jalapic/6ca3ece44bdcdc522bb735f183aa0ca0/raw/1a07f469eff08117121b6cbebcd22fb7569e3ee8/compete_extra.R")
source('lib/matrixPlotParameter.R')
source('lib/residualStructureScore.R')





#### then load in the data here ###

full_data <- readRDS("data/mice.RData")
# A=c9, B=c10, C=c12, D=c15, E=c16, F=c17, G=c18, H=c37, I=c38. J=c45
cohort_names <- paste("cohort",c(9,10,12,15,16,17,18,37,38,45),sep='')
cohort_short_names <- paste("C",c(9,10,12,15,16,17,18,37,38,45),sep='')
cut_off <- 3
mice_number <- 12


model1_fn <- list(alpha.fun = function(x,y,eta1,eta2,eta3){return(eta1*x*y*exp(-eta2*abs(x-y))/(1+exp(-eta3*(x-y))))})

model3_fn <- list(alpha.fun = function(x,y,eta1,eta2){return(eta1*x*y*exp(-eta2*abs(x-y)))},
                  q1.fun = function(x,y,eta3){return(exp(-eta3*x))},
                  q0.fun = function(x,y,eta3){return(exp(-eta3*y))})

# Define the cohorts will be fitted
fit_cohorts <- c(1:10)
naive_rank_10 <- list()
expert_rank_10 <- list()
for(current_cohort in fit_cohorts){
  naive_rank_10[[current_cohort]] <- naiveRankHierarchy(full_data[[cohort_names[current_cohort]]])
  expert_rank_10[[current_cohort]] <- expertRankHierarchy(full_data[[cohort_names[current_cohort]]])
}




current_cohort <- cohort_id
print(paste("Cohort",current_cohort))


# then fit I-MMHP here for each cohort ####

print(current_cohort)
stan_input_lst <- prepareDataStanIMMHP(current_cohort)

fit_mmhp_sep <- stan("lib/I-MMHP.stan",
                     data=stan_input_lst,
                     warmup = 2000, iter = 3000, thin = 4, chains = 4) 
sim_mmhp_sep <- rstan::extract(fit_mmhp_sep)
dir.create(paste(save_data_path, cohort_names[current_cohort],sep=''), recursive = TRUE, showWarnings = FALSE)
save(sim_mmhp_sep, fit_mmhp_sep,
     file = paste(save_data_path,cohort_names[current_cohort],
                  "/sep_mmhp_stan_result_",cohort_names[current_cohort],
                  ".RData",sep=''))

# Interpolate Latent States ####

print(current_cohort)
load(paste(save_data_path,cohort_names[current_cohort],
           "/sep_mmhp_stan_result_",cohort_names[current_cohort],
           ".RData",sep=''))
clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
return_df <- cleanObservationPeriod(current_cohort, clean_data)
unique_pairs_df <- return_df %>% group_by(initiator, recipient) %>%
  summarize(count=n(),
            observe=list(observe.id),
            observe.length=list(observe.time),
            no.events=list(no.events))
unique_observe_win <- unique(return_df[,c("observe.id","observe.time")])

num_winds <- nrow(unique_observe_win)

state_array_list <- list() 
initial_state_list <- list()
termination_state_list <- list()
interpolation_array_list <- list()
#no_segments <- 5000

for(pair in 1:nrow(unique_pairs_df)){
  print(pair)
  current_initiator <- as.numeric(unique_pairs_df[pair,"initiator"])
  current_recipient <- as.numeric(unique_pairs_df[pair,"recipient"])
  current_window_vec <- unique_pairs_df$observe[[pair]]
  
  param <- rep(list(),1000)
  delta_vec <- rep(0,1000)
  for(current_sim in 1:1000){
    param[[current_sim]] <- list(lambda1=sim_mmhp_sep$lambda1[current_sim,pair],
                                 lambda0=sim_mmhp_sep$lambda0[current_sim,pair],
                                 alpha=sim_mmhp_sep$alpha[current_sim,pair],
                                 beta=sim_mmhp_sep$beta[current_sim,pair],
                                 q1=sim_mmhp_sep$q1[current_sim,pair],
                                 q2=sim_mmhp_sep$q2[current_sim,pair])
    delta_vec[current_sim] <- sim_mmhp_sep$delta_1[current_sim,pair]
  }
  
  state_array_list[[pair]] <- list()
  initial_state_list[[pair]] <- list()
  termination_state_list[[pair]] <- list()
  interpolation_array_list[[pair]] <- list()
  
  for(current_win in 1:num_winds){
    
    row_indicator <- return_df$initiator==current_initiator&return_df$recipient==current_recipient&return_df$observe.id==current_win
    
    if(current_win %in% current_window_vec){
      time_vec <- return_df[row_indicator,"event.times"][[1]]
      observe_period <- unique_observe_win[unique_observe_win$observe.id==current_win,"observe.time"]
      state_array_list[[pair]][[current_win]] <- matrix(0,nrow=length(time_vec),ncol=1000)
    }
    else {
      time_vec <- NULL
      observe_period <- return_df[return_df$observe.id == current_win,"observe.time"][1]
      state_array_list[[pair]][[current_win]] <- matrix(0,nrow=1,ncol=1000)
      # is this the right length for state_array_list?
    }
    
    
    #time_vec <- return_df[row_indicator,"event.times"][[1]]
    #observe_period <- unique_observe_win[unique_observe_win$observe.id==current_win,"observe.time"]
    time_segment <- seq(0,observe_period,length.out=no_segments)
    
    #state_array_list[[pair]][[current_win]] <- matrix(0,nrow=length(time_vec),ncol=1000)
    initial_state_list[[pair]][[current_win]] <- matrix(0,nrow=1,ncol=1000)
    termination_state_list[[pair]][[current_win]] <- matrix(0,nrow=1,ncol=1000)
    interpolation_array_list[[pair]][[current_win]] <- matrix(0,nrow=no_segments,ncol=1000)
    
    for(current_sim in 1:1000){
      ## latent states at event times
      viterbi_result <- myViterbiWithInitial(events = time_vec, param = param[[current_sim]],
                                             initial.p = delta_vec[current_sim],
                                             termination = observe_period)
      state_array_list[[pair]][[current_win]][,current_sim] <- viterbi_result$zt_v
      initial_state_list[[pair]][[current_win]][1,current_sim] <- viterbi_result$initial_state
      termination_state_list[[pair]][[current_win]][1,current_sim] <- viterbi_result$termination_state
      
      ## interpolation
      latent_inter <- interpolateLatentTrajectory(param[[current_sim]], time_vec, viterbi_result$zt_v,
                                                  initial.state = viterbi_result$initial_state,
                                                  termination.time=observe_period,
                                                  termination.state = viterbi_result$termination_state)
      step_fun_est <- stepfun(latent_inter$x.hat,2-latent_inter$z.hat)
      interpolation_array_list[[pair]][[current_win]][,current_sim] <- step_fun_est(time_segment)
    }
  }
}
save(state_array_list,initial_state_list,termination_state_list,
     interpolation_array_list,no_segments,
     file=paste(save_data_path,cohort_names[current_cohort],
                "/mmhp_est_zt_",cohort_names[current_cohort],".RData",sep=''))

#### predictions ####

print(current_cohort)
stan_train_input_lst <- prepareDataStanTrainIMMHP(current_cohort)

fit_mmhp_sep <- stan("lib/I-MMHP.stan",
                     data=stan_train_input_lst,
                     warmup = 1000, iter = 2000, thin = 4, chains = 4) 
sim_mmhp_sep <- rstan::extract(fit_mmhp_sep)
dir.create(paste(save_data_path, cohort_names[current_cohort],sep=''), 
           recursive = TRUE, showWarnings = FALSE)
save(fit_mmhp_sep, sim_mmhp_sep,
     file = paste(save_data_path,cohort_names[current_cohort],
                  "/predict_immhp_stan_result_",cohort_names[current_cohort],
                  ".RData",sep=''))


# then the pearson residuals for this fit ####
mice_number <- 12

clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
return_df <- cleanObservationPeriod(current_cohort, clean_data)
unique_pairs_df <- return_df %>% group_by(initiator, recipient) %>%
  summarize(count=n(),
            observe=list(observe.id),
            observe.length=list(observe.time),
            no.events=list(no.events))
unique_observe_win <- unique(return_df[,c("observe.id","observe.time")])

num_winds <- nrow(unique_observe_win)


load(paste(save_data_path, cohort_names[current_cohort], 
           "/sep_mmhp_stan_result_", cohort_names[current_cohort],".RData",sep=''))
load(paste(save_data_path,cohort_names[current_cohort],
           "/mmhp_est_zt_",cohort_names[current_cohort],".RData",sep=''))

mmhp_residual_matrix <- matrix(0,ncol=mice_number,nrow=mice_number)

mmhp_residual_array <- array(0,dim = c(mice_number,mice_number,num_winds))

for(i in 1:mice_number){
  print(i)
  for(j in 1:mice_number){
    pair <- which(unique_pairs_df$initiator==i&unique_pairs_df$recipient==j)
    if(length(pair)>0&(i!=j)){
      par_est <- list(lambda0=mean(sim_mmhp_sep$lambda0[,pair]),
                      lambda1=mean(sim_mmhp_sep$lambda1[,pair]),
                      alpha=mean(sim_mmhp_sep$alpha[,pair]),
                      beta=mean(sim_mmhp_sep$beta[,pair]),
                      q1=mean(sim_mmhp_sep$q1[,pair]),
                      q2=mean(sim_mmhp_sep$q2[,pair]))
      current_window_vec <- unique_pairs_df$observe[[pair]]
      all_residual <- 0
      for(cur in c(1:num_winds)){ ## check length > 2
        if(cur %in% current_window_vec) {
          cur_win <- cur#current_window_vec[cur]
          current_event_time <- return_df[return_df$initiator==i&
                                            return_df$recipient==j&
                                            return_df$observe.id==cur_win,"event.times"][[1]]
          # I think this just returns the windows where there are events
          
          current_obs_time <- return_df[return_df$initiator==i&
                                          return_df$recipient==j&
                                          return_df$observe.id==cur_win,"observe.time"]
          time_segment <- seq(0,current_obs_time,length.out=no_segments)
          latent_mean <- apply(interpolation_array_list[[pair]][[cur_win]],1,mean)
          latent_event <- as.numeric(apply(2-state_array_list[[pair]][[cur_win]],1,mean) > 0.5) 
          est.intensity <- mmhpIntensityNumeric(params=par_est,
                                                t=current_event_time,
                                                time.vec=time_segment,
                                                latent.vec=latent_mean)
          est.intensity.events <- mmhpIntensityAtEvents(params=par_est, t=current_event_time,
                                                        latent_z=latent_event)
          all_residual <- all_residual + sum(1/sqrt(est.intensity.events))-
            sum(sqrt(est.intensity))*(time_segment[2]-time_segment[1])
          mmhp_residual_array[i,j,cur] <- sum(1/sqrt(est.intensity.events))-
            sum(sqrt(est.intensity))*(time_segment[2]-time_segment[1])
          
        }
        else {
          current_event_time <- NULL
          # then need to get the current_obs_time
          current_obs_time <- return_df[return_df$observe.id == cur,"observe.time"][1]
          time_segment <- seq(0,current_obs_time,length.out=no_segments)
          latent_mean <- apply(interpolation_array_list[[pair]][[cur]],1,mean)
          latent_event <- as.numeric(apply(2-state_array_list[[pair]][[cur]],1,mean) > 0.5) 
          est.intensity <- mmhpIntensityNumeric_win(params=par_est,
                                                    t=current_event_time,
                                                    time.vec=time_segment,
                                                    latent.vec=latent_mean)
          # est.intensity.events <- mmhpIntensityAtEvents(params=par_est, t=current_event_time,
          #                                               latent_z=latent_event)
          all_residual <- all_residual -
            sum(sqrt(est.intensity))*(time_segment[2]-time_segment[1])
          mmhp_residual_array[i,j,cur] <- -1*sum(sqrt(est.intensity))*(time_segment[2]-time_segment[1])
        }
        # cur_win <- current_window_vec[cur]
        # current_event_time <- return_df[return_df$initiator==i&
        #                                   return_df$recipient==j&
        #                                   return_df$observe.id==cur_win,"event.times"][[1]]
        # current_obs_time <- return_df[return_df$initiator==i&
        #                                 return_df$recipient==j&
        #                                 return_df$observe.id==cur_win,"observe.time"]
        # time_segment <- seq(0,current_obs_time,length.out=no_segments)
        # latent_mean <- apply(interpolation_array_list[[pair]][[cur_win]],1,mean)
        # latent_event <- as.numeric(apply(2-state_array_list[[pair]][[cur_win]],1,mean) > 0.5) 
        # est.intensity <- mmhpIntensityNumeric(params=par_est,
        #                                       t=current_event_time,
        #                                       time.vec=time_segment,
        #                                       latent.vec=latent_mean)
        # est.intensity.events <- mmhpIntensityAtEvents(params=par_est, t=current_event_time,
        #                                               latent_z=latent_event)
        # all_residual <- all_residual + sum(1/sqrt(est.intensity.events))-
        #   sum(sqrt(est.intensity))*(time_segment[2]-time_segment[1])
        # 
        # mmhp_residual_array[i,j,cur]
      }
      mmhp_residual_matrix[i,j] <- all_residual
    }
  }
}


# then save this mmhp_residual_matrix
saveRDS(mmhp_residual_matrix,
        file = paste(save_data_path,cohort_names[current_cohort],
                     "/immhp_pr_matrix_",cohort_names[current_cohort],
                     ".RDS",sep=''))


