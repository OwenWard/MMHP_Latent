### Diagnostics at a cohort level ####

### Construct cohort level diagnostic statistics for model fit for each of the
### 4 models (c-hp,dc-hp,c-mmhp,i-mmhp) fit in this paper


### code ###
## run this if running on the cluster
source("/rigel/stats/users/ogw2103/code/MMHP/MMHP_Latent/run_scripts/cluster_setup.R")
### set cohort_id based on job num
jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
cohort_id <- jobid
# cohort_id <- 1
####
data_path <- "output/"


# library(rstan)
# options(mc.cores = parallel::detectCores())


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






#### then load in the data ####

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

clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
return_df <- cleanObservationPeriod(current_cohort, clean_data)
unique_pairs_df <- return_df %>% group_by(initiator, recipient) %>%
  summarize(count=n(),
            observe=list(observe.id),
            observe.length=list(observe.time),
            no.events=list(no.events))
unique_observe_win <- unique(return_df[,c("observe.id","observe.time")])

num_winds <- nrow(unique_observe_win)



#### Diagnostics M1 Cohort Hawkes ####
raw_real_c_hawkes_vec  <- numeric(0)
pr_real_c_hawkes_vec <- numeric(0)
Lambda_c_hawkes_matrix <- matrix(rep(list(),mice_number*mice_number),
                                 nrow=mice_number, ncol=mice_number)


load(paste(data_path, cohort_names[current_cohort],
           "/cohort_hp_stan_result_", cohort_names[current_cohort],".RData",sep=''))

# create a list to store the hawkes parameters for this
model1_par_est <- list(lambda0=mean(sim_cohort_hp$lambda0),
                       eta_1=mean(sim_cohort_hp$eta_1),
                       eta_2=mean(sim_cohort_hp$eta_2),
                       eta_3=mean(sim_cohort_hp$eta_3),
                       beta=mean(sim_cohort_hp$beta),
                       f=apply(sim_cohort_hp$f,2,mean))
model1_par_matrix <- list(lambda0_matrix=matrix(model1_par_est$lambda0,
                                                nrow=mice_number,ncol=mice_number),
                          alpha_matrix=formMatrix(function(x,y) model1_fn$alpha.fun(x,y,
                                                                                    model1_par_est$eta_1,
                                                                                    model1_par_est$eta_2,
                                                                                    model1_par_est$eta_3),
                                                  model1_par_est$f),
                          beta_matrix=matrix(model1_par_est$beta,
                                             nrow=mice_number,ncol=mice_number))



for(i in 1:mice_number){
  print(i)
  for(j in 1:mice_number){
    pair <- which(unique_pairs_df$initiator==i&unique_pairs_df$recipient==j)
    if(length(pair)>0&(i!=j)){
      current_window_vec <- unique_pairs_df$observe[[pair]]
      all_rescaled_interevent <- numeric(0)
      for(cur in c(1:length(current_window_vec))){ 
        cur_win <- current_window_vec[cur]
        current_event_time <- return_df[return_df$initiator==i&
                                          return_df$recipient==j&
                                          return_df$observe.id==cur_win,"event.times"][[1]]
        observe_period <- return_df[return_df$initiator==i&
                                      return_df$recipient==j&
                                      return_df$observe.id==cur_win,"observe.time"]
        #### update this from the stan fit
        hawkes.par <- c(model1_par_matrix$lambda0_matrix[i,j],model1_par_matrix$alpha[i,j],
                        model1_par_matrix$beta[i,j])
        # hawkes.par <- optim(par=c(1,0,100), fn=uniHawkesNegLogLik, method="CG", 
        #                     t=current_event_time, termination = observe_period)
        all_Lambda <- uniHawkesIntegralIntensity(object=list(lambda0=hawkes.par[1],
                                                             alpha=hawkes.par[2],
                                                             beta=hawkes.par[3]),
                                                 events=current_event_time,
                                                 termination = observe_period)
        all_prresidual <- uniHawkesPearsonResidual(object=list(lambda0=hawkes.par[1],
                                                               alpha=abs(hawkes.par[2]),
                                                               beta=hawkes.par[3]),
                                                   events=current_event_time,
                                                   termination = observe_period)
        
        ## rescaled interevent time
        hawkes_obj <- list(lambda0 = hawkes.par[1],
                           alpha = hawkes.par[2],
                           beta = hawkes.par[3])
        # Lambda_vec <- uniHawkesCompensator(lambda0=hawkes.par[1],
        #                                    alpha=hawkes.par[2],
        #                                    beta=hawkes.par[3],
        #                                    c(0,current_event_time))
        Lambda_vec <- uniHawkesCompensator(hawkes_obj,current_event_time)
        ## should this be c(0,events)?
        
        all_rescaled_interevent <- c(all_rescaled_interevent,Lambda_vec)
        raw_real_c_hawkes_vec <- c(raw_real_c_hawkes_vec,all_Lambda)
        pr_real_c_hawkes_vec <- c(pr_real_c_hawkes_vec,all_prresidual)
      }
      Lambda_c_hawkes_matrix[i,j][[1]] <- all_rescaled_interevent
    }
  }
}




#### Diagnostics M2 Cohort Degree Corrected Hawkes ####
raw_real_dc_hawkes_vec  <- numeric(0)
pr_real_dc_hawkes_vec <- numeric(0)
Lambda_dc_hawkes_matrix <- matrix(rep(list(),mice_number*mice_number),
                                  nrow=mice_number, ncol=mice_number)


load(paste(data_path, cohort_names[current_cohort],
           "/cohort_dchp_stan_result_", cohort_names[current_cohort],".RData",sep=''))

# create a list to store the hawkes parameters for this
model2_par_est <- list(gamma = apply(sim_cohort_dchp$gamma,2,mean),
                       zeta = apply(sim_cohort_dchp$zeta,2,mean),
                       #lambda0=mean(sim_cohort_dchp$lambda0),
                       eta_1=mean(sim_cohort_dchp$eta_1),
                       eta_2=mean(sim_cohort_dchp$eta_2),
                       eta_3=mean(sim_cohort_dchp$eta_3),
                       beta=mean(sim_cohort_dchp$beta),
                       f=apply(sim_cohort_dchp$f,2,mean))
model2_par_matrix <- list(lambda0_matrix= outer(model2_par_est$gamma,
                                                model2_par_est$zeta,FUN = "+"),
                          alpha_matrix=formMatrix(function(x,y) model1_fn$alpha.fun(x,y,
                                                                                    model2_par_est$eta_1,
                                                                                    model2_par_est$eta_2,
                                                                                    model2_par_est$eta_3),
                                                  model2_par_est$f),
                          beta_matrix=matrix(model2_par_est$beta,
                                             nrow=mice_number,ncol=mice_number))



for(i in 1:mice_number){
  print(i)
  for(j in 1:mice_number){
    pair <- which(unique_pairs_df$initiator==i&unique_pairs_df$recipient==j)
    if(length(pair)>0&(i!=j)){
      current_window_vec <- unique_pairs_df$observe[[pair]]
      all_rescaled_interevent <- numeric(0)
      for(cur in c(1:length(current_window_vec))){ 
        cur_win <- current_window_vec[cur]
        current_event_time <- return_df[return_df$initiator==i&
                                          return_df$recipient==j&
                                          return_df$observe.id==cur_win,"event.times"][[1]]
        observe_period <- return_df[return_df$initiator==i&
                                      return_df$recipient==j&
                                      return_df$observe.id==cur_win,"observe.time"]
        #### update this from the stan fit
        hawkes.par <- c(model2_par_matrix$lambda0_matrix[i,j],model2_par_matrix$alpha[i,j],
                        model2_par_matrix$beta[i,j])
        # hawkes.par <- optim(par=c(1,0,100), fn=uniHawkesNegLogLik, method="CG", 
        #                     t=current_event_time, termination = observe_period)
        all_Lambda <- uniHawkesIntegralIntensity(object=list(lambda0=hawkes.par[1],
                                                             alpha=hawkes.par[2],
                                                             beta=hawkes.par[3]),
                                                 events=current_event_time,
                                                 termination = observe_period)
        all_prresidual <- uniHawkesPearsonResidual(object=list(lambda0=hawkes.par[1],
                                                               alpha=abs(hawkes.par[2]),
                                                               beta=hawkes.par[3]),
                                                   events=current_event_time,
                                                   termination = observe_period)
        
        ## rescaled interevent time
        hawkes_obj <- list(lambda0 = hawkes.par[1],
                           alpha = hawkes.par[2],
                           beta = hawkes.par[3])
        Lambda_vec <- uniHawkesCompensator(hawkes_obj,
                                           current_event_time)
        all_rescaled_interevent <- c(all_rescaled_interevent,Lambda_vec)
        raw_real_dc_hawkes_vec <- c(raw_real_dc_hawkes_vec,all_Lambda)
        pr_real_dc_hawkes_vec <- c(pr_real_dc_hawkes_vec,all_prresidual)
      }
      Lambda_dc_hawkes_matrix[i,j][[1]] <- all_rescaled_interevent
    }
  }
}


## Diagnostics M3 CMMHP ####
load(paste(data_path, cohort_names[current_cohort],
           "/cohort_mmhp_stan_result_", cohort_names[current_cohort],".RData",sep=''))
load(paste(data_path,cohort_names[current_cohort],
           "/cmmhp_est_zt_",cohort_names[current_cohort],".RData",sep=''))


raw_real_cmmhp_vec  <- numeric(0)
pr_real_cmmhp_vec <- numeric(0)
Lambda_cmmhp_matrix <- matrix(rep(list(),mice_number*mice_number),
                              nrow=mice_number, ncol=mice_number)

real_N_vec <- numeric(0)
for(i in 1:mice_number){
  print(i)
  for(j in 1:mice_number){
    pair <- which(unique_pairs_df$initiator==i&unique_pairs_df$recipient==j)
    if(length(pair)>0&(i!=j)){
      par_est <- list(lambda0=mean(sim_cohort_mmhp$lambda0[,pair]),
                      lambda1=mean(sim_cohort_mmhp$lambda1[,pair]),
                      alpha=mean(sim_cohort_mmhp$alpha[,pair]),
                      beta=mean(sim_cohort_mmhp$beta),
                      q1=mean(sim_cohort_mmhp$q1[,pair]),
                      q2=mean(sim_cohort_mmhp$q2[,pair]))
      current_window_vec <- unique_pairs_df$observe[[pair]]
      all_rescaled_interevent <- numeric(0)
      #for(cur in c(1:length(current_window_vec))){ 
      for(cur in current_window_vec){
        #cur_win <- current_window_vec[cur]
        cur_win <- cur
        current_event_time <- return_df[return_df$initiator==i&
                                          return_df$recipient==j&
                                          return_df$observe.id==cur_win,"event.times"][[1]]
        current_obs_time <- return_df[return_df$initiator==i&
                                        return_df$recipient==j&
                                        return_df$observe.id==cur_win,"observe.time"]
        time_segment <- seq(0,current_obs_time,length.out=no_segments)
        ### I need to adjust the indexing of these to get the corresponding windows,
        ### below was written for only windows with events, when state array for all now
        latent_mean <- apply(interpolation_array_list[[pair]][[cur_win]],1,mean)
        latent_event <- as.numeric(apply(2-state_array_list[[pair]][[cur_win]],1,mean) > 0.5) 
        
        ## Pearson
        est.intensity <- mmhpIntensityNumeric(params=par_est,
                                              t=current_event_time,
                                              time.vec=time_segment,
                                              latent.vec=latent_mean)
        # this is causing a bug...
        
        est.intensity.events <- mmhpIntensityAtEvents(params=par_est, t=current_event_time,
                                                      latent_z=latent_event)
        all_prresidual <- sum(1/sqrt(est.intensity.events))-
          sum(sqrt(est.intensity))*(time_segment[2]-time_segment[1])
        all_Lambda <- sum(est.intensity)*(time_segment[2]-time_segment[1])
        
        ## rescaled interevent time
        Lambda_vec <- rep(0,length(current_event_time)+1)
        temp_time_vec <- c(0,current_event_time,current_obs_time)
        for(m in 1:(length(current_event_time)+1)){
          interevent_idx <- (time_segment>=temp_time_vec[m]) & (time_segment<temp_time_vec[m+1])
          Lambda_vec[m] <- sum(est.intensity[interevent_idx])*(time_segment[2]-time_segment[1])
        }
        real_N_vec <- c(real_N_vec, length(current_event_time))
        raw_real_cmmhp_vec <- c(raw_real_cmmhp_vec,all_Lambda)
        pr_real_cmmhp_vec <- c(pr_real_cmmhp_vec,all_prresidual)
        all_rescaled_interevent <- c(all_rescaled_interevent,Lambda_vec)
      }  
      Lambda_cmmhp_matrix[i,j][[1]] <- all_rescaled_interevent
    }
  }
}

rm(interpolation_array_list)
rm(state_array_list)
rm(termination_state_list)
rm(initial_state_list)


#### Diagnostics I-MMHP ####
load(paste(data_path, cohort_names[current_cohort],
           "/sep_mmhp_stan_result_", cohort_names[current_cohort],".RData",sep=''))
load(paste(data_path,cohort_names[current_cohort],
           "/mmhp_est_zt_",cohort_names[current_cohort],".RData",sep=''))

## indep - mmhp
raw_real_mmhp_vec  <- numeric(0)
pr_real_mmhp_vec <- numeric(0)
Lambda_mmhp_matrix <- matrix(rep(list(),mice_number*mice_number),
                             nrow=mice_number, ncol=mice_number)

real_N_vec <- numeric(0)
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
      all_rescaled_interevent <- numeric(0)
      for(cur in c(1:length(current_window_vec))){ 
        cur_win <- current_window_vec[cur]
        current_event_time <- return_df[return_df$initiator==i&
                                          return_df$recipient==j&
                                          return_df$observe.id==cur_win,"event.times"][[1]]
        current_obs_time <- return_df[return_df$initiator==i&
                                        return_df$recipient==j&
                                        return_df$observe.id==cur_win,"observe.time"]
        time_segment <- seq(0,current_obs_time,length.out=no_segments)
        latent_mean <- apply(interpolation_array_list[[pair]][[cur_win]],1,mean)
        latent_event <- as.numeric(apply(2-state_array_list[[pair]][[cur_win]],1,mean) > 0.5) 
        
        ## Pearson
        est.intensity <- mmhpIntensityNumeric(params=par_est,
                                              t=current_event_time,
                                              time.vec=time_segment,
                                              latent.vec=latent_mean)
        est.intensity.events <- mmhpIntensityAtEvents(params=par_est, t=current_event_time,
                                                      latent_z=latent_event)
        all_prresidual <- sum(1/sqrt(est.intensity.events))-
          sum(sqrt(est.intensity))*(time_segment[2]-time_segment[1])
        all_Lambda <- sum(est.intensity)*(time_segment[2]-time_segment[1])
        
        ## rescaled interevent time
        Lambda_vec <- rep(0,length(current_event_time)+1)
        temp_time_vec <- c(0,current_event_time,current_obs_time)
        for(m in 1:(length(current_event_time)+1)){
          interevent_idx <- (time_segment>=temp_time_vec[m]) & (time_segment<temp_time_vec[m+1])
          Lambda_vec[m] <- sum(est.intensity[interevent_idx])*(time_segment[2]-time_segment[1])
        }
        real_N_vec <- c(real_N_vec, length(current_event_time))
        raw_real_mmhp_vec <- c(raw_real_mmhp_vec,all_Lambda)
        pr_real_mmhp_vec <- c(pr_real_mmhp_vec,all_prresidual)
        all_rescaled_interevent <- c(all_rescaled_interevent,Lambda_vec)
      }  
      Lambda_mmhp_matrix[i,j][[1]] <- all_rescaled_interevent
    }
  }
}

#### Save Output ####
### save all these outputs in a nice format
save(raw_real_c_hawkes_vec, pr_real_c_hawkes_vec,Lambda_c_hawkes_matrix,
     raw_real_dc_hawkes_vec, pr_real_dc_hawkes_vec, Lambda_dc_hawkes_matrix,
     raw_real_cmmhp_vec, pr_real_cmmhp_vec, Lambda_cmmhp_matrix,
     raw_real_mmhp_vec, pr_real_mmhp_vec,Lambda_mmhp_matrix,real_N_vec,
     file = paste(data_path, cohort_names[current_cohort],
                  "/real_diag_four_models_",
                  cohort_names[current_cohort],".RData",sep='') )


