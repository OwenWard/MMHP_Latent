#### Fit the dsnl model ####


### code ###
## run this if running on the cluster
source("/rigel/stats/users/ogw2103/code/MMHP/MMHP_Latent/run_scripts/cluster_setup.R")
# ### set cohort_id based on job num
jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
cohort_id <- jobid
# cohort_id <- 1
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






#### then load in the data ###

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


#### Fit DSNL stan Model ####

current_cohort <- cohort_id
print(paste("Cohort",current_cohort))

# prepare input data for stan
no_days <- length(unique(full_data[[cohort_names[current_cohort]]]$day))
N_count_day <- array(0,dim=c(no_days,12,12))
N_count_indicator_day <- array(0,dim=c(no_days,12,12))

for(k in c(1:no_days)){
  temp_matrix <- matrix(0,12,12)
  if(k==1){
    N_count_day[k,,] <- cleanDataForTraining(all_raw_data = full_data[[cohort_names[current_cohort]]], 
                                             cut_off = 0, N = 12, train_period = 1)$N_count
  }else{
    clean_data_test_1 <- cleanDataForTraining(all_raw_data = full_data[[cohort_names[current_cohort]]], 
                                              cut_off = 0, N = 12, train_period = k-1)
    clean_data_test_2 <- cleanDataForTraining(all_raw_data = full_data[[cohort_names[current_cohort]]], 
                                              cut_off = 0, N = 12, train_period = k)
    N_count_day[k,,] <- (clean_data_test_2$N_count - clean_data_test_1$N_count)
    temp_matrix[N_count_day[k,,]>0] <- 1
    N_count_indicator_day[k,,] <- temp_matrix
  }
}

# fit in stan
fit_dsnl <- stan("lib/dsnl_poisson.stan",
                 data=list(day=no_days,
                           Gt=N_count_day,
                           c=1,
                           sigma=1),
                 iter=1000, chains=4, control=list(adapt_delta=0.99))
sim_dsnl <- rstan::extract(fit_dsnl)
save(sim_dsnl, fit_dsnl,
     file = paste(save_data_path,cohort_names[current_cohort],
                  "/dsnl_poisson_stan_result_",cohort_names[current_cohort],
                  ".RData",sep=''))

#### Predict DSNL Model ####

no_days <- length(unique(full_data[[cohort_names[current_cohort]]]$day))
no_days = 15
N_count_day <- array(0,dim=c(no_days,12,12))
N_count_indicator_day <- array(0,dim=c(no_days,12,12))

for(k in c(1:no_days)){
  temp_matrix <- matrix(0,12,12)
  if(k==1){
    N_count_day[k,,] <- cleanDataForTraining(all_raw_data = full_data[[cohort_names[current_cohort]]], 
                                             cut_off = 0, N = 12, train_period = 1)$N_count
  }else{
    clean_data_test_1 <- cleanDataForTraining(all_raw_data = full_data[[cohort_names[current_cohort]]], 
                                              cut_off = 0, N = 12, train_period = k-1)
    clean_data_test_2 <- cleanDataForTraining(all_raw_data = full_data[[cohort_names[current_cohort]]], 
                                              cut_off = 0, N = 12, train_period = k)
    N_count_day[k,,] <- (clean_data_test_2$N_count - clean_data_test_1$N_count)
    temp_matrix[N_count_day[k,,]>0] <- 1
    N_count_indicator_day[k,,] <- temp_matrix
  }
}


no_days <- length(unique(full_data[[cohort_names[current_cohort]]]$day))
N_count_train <- cleanDataForTraining(all_raw_data = full_data[[cohort_names[current_cohort]]], 
                                      cut_off = cut_off, N = 12, train_period = 15)$N_count

fit_predict_dsnl <- stan("lib/dsnl_predict.stan",
                         data=list(day=15,
                                   Gt=N_count_day,
                                   c=1,
                                   sigma=1),
                         iter=1500, chains=4, control=list(adapt_delta=0.999))



sim_predict_dsnl <- rstan::extract(fit_predict_dsnl)
dir.create(paste(save_data_path, cohort_names[current_cohort],sep=''), 
           recursive = TRUE, showWarnings = FALSE)
save(fit_predict_dsnl, sim_predict_dsnl,
     file = paste(save_data_path,cohort_names[current_cohort],
                  "/predict_dsnl_poisson_stan_result_",cohort_names[current_cohort],
                  ".RData",sep=''))
