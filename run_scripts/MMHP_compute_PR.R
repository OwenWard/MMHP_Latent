# file to work on estimating the Pearson residuals 
# for windows with no events

# this file refits the complete C-MMHP model
# including prediction and Pearson residuals for each cohort
#

### code ###
## run this if running on the cluster
#source("/rigel/stats/users/ogw2103/code/MMHP/MMHP_Latent/run_scripts/cluster_setup.R")
# ### set cohort_id based on job num
# jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# jobid <- as.numeric(jobid)
# cohort_id <- jobid
cohort_id <- 1
####
save_data_path <- "output/"

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




#Jing_data_path <- "/rigel/stats/users/ogw2103/code/MMHP/Latent_Sims/real_data/"
save_data_path <- "/rigel/stats/users/ogw2103/code/MMHP/Latent_Sims/paper_results/"
# to run locally
#save_data_path <- "C:/Users/owenw/Google Drive/Tian/Research+JingWu/thesis_data/part2/real_data/"

#### then load in the stan fit here ###

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


## testing locally
cohort_id <- 1
####

current_cohort <- cohort_id
print(paste("Cohort",current_cohort))


# then fit C-MMHP here for each cohort
#######
#######
#######
## The stan model used here is important
#######
#######
#######


print(paste("Cohort",current_cohort))

# then the pearson residuals for this fit will also be computed here
mice_number <- 12

load(paste(save_data_path,cohort_names[current_cohort],
           "/cohort_mmhp_stan_result_",cohort_names[current_cohort],
           ".RData",sep=''))
load(paste(save_data_path,cohort_names[current_cohort],
           "/cmmhp_est_zt_",cohort_names[current_cohort],".RData",sep=''))

clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
return_df <- cleanObservationPeriod(current_cohort, clean_data)
unique_pairs_df <- return_df %>% group_by(initiator, recipient) %>%
  summarize(count=n(),
            observe=list(observe.id),
            observe.length=list(observe.time),
            no.events=list(no.events))
unique_observe_win <- unique(return_df[,c("observe.id","observe.time")])

num_windows <- nrow(unique_observe_win) # this is safe?


model3_par_est <- list(lambda0=mean(sim_cohort_mmhp$lambda0),
                       lambda1=mean(sim_cohort_mmhp$lambda1),
                       eta_1=mean(sim_cohort_mmhp$eta_1),
                       eta_2=mean(sim_cohort_mmhp$eta_2),
                       eta_3=mean(sim_cohort_mmhp$eta_3),
                       beta=mean(sim_cohort_mmhp$beta),
                       f=apply(sim_cohort_mmhp$f,2,mean))
model3_par_matrix <- list(lambda0_matrix=matrix(model3_par_est$lambda0,
                                                nrow=mice_number,ncol=mice_number),
                          lambda1_matrix=matrix(model3_par_est$lambda1,
                                                nrow=mice_number,ncol=mice_number),
                          alpha_matrix=formMatrix(function(x,y) model3_fn$alpha.fun(x,y,
                                                                                    model3_par_est$eta_1,
                                                                                    model3_par_est$eta_2),
                                                  model3_par_est$f),
                          beta_matrix=matrix(model3_par_est$beta,
                                             nrow=mice_number,ncol=mice_number),
                          q1_matrix=formMatrix(function(x,y) model3_fn$q1.fun(x,y,model3_par_est$eta_3),
                                               model3_par_est$f),
                          q2_matrix=formMatrix(function(x,y) model3_fn$q0.fun(x,y,model3_par_est$eta_3),
                                               model3_par_est$f))
m3_residual_matrix <- matrix(0,ncol=mice_number,nrow=mice_number)
no_segments <- 500 # changed from 5000
window_pr <- c()
for(i in 1:mice_number){
  print(i)
  for(j in 1:mice_number){
    pair <- which(unique_pairs_df$initiator==i&unique_pairs_df$recipient==j)
    if(length(pair)>0&(i!=j)){
      par_est <- lapply(model3_par_matrix, function(x) x[i,j])
      names(par_est) <- c("lambda0","lambda1","alpha","beta","q1","q2")
      current_window_vec <- unique_pairs_df$observe[[pair]]
      all_residual <- 0
      for(cur in c(1:num_windows)) { ## check length > 2
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
          
        }
        
        # print(sum(1/sqrt(est.intensity.events))-
        #         sum(sqrt(est.intensity))*(time_segment[2]-time_segment[1]))
        # print(summary(latent_mean))
        # window_pr <- c(window_pr,sum(1/sqrt(est.intensity.events))-
        #                  sum(sqrt(est.intensity))*(time_segment[2]-time_segment[1]))
      }
      m3_residual_matrix[i,j] <- all_residual
    }
  }
}


cat("Finished without errors")
# then save this mmhp_residual_matrix
# saveRDS(m3_residual_matrix,
#         file = paste(save_data_path,cohort_names[current_cohort],
#                      "/cmmhp_pr_matrix_",cohort_names[current_cohort],
#                      ".RDS",sep=''))


