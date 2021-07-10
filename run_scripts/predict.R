### Predictions at a cohort level ####


### code ###
## run this if running on the cluster
source("/moto/stats/users/ogw2103/Code/MMHP_Latent/run_scripts/cluster_setup.R")
### set cohort_id based on job num
jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
cohort_id <- jobid
#cohort_id <- 1
####
data_path <- "output/revisions/"


library(rstan)
options(mc.cores = parallel::detectCores())


library(compete)
#library(RColorBrewer)
#library(Hmisc)
#library(wCorr)
#library(tidyverse)
library(dplyr)
library(R.utils)
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

predict_sim <- 1000

print(current_cohort)
clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
return_df <- cleanObservationPeriod(current_cohort, clean_data)
to_predice_obs <- unique(return_df[(return_df$day>=16)&(return_df$day<=21),c("observe.id","observe.time")])

train_fit <- prepareDataStanTrain(current_cohort)


#### Predict M1 ####
load(paste(data_path,cohort_names[current_cohort],
           "/cohort_hp_predict_stan_result_",cohort_names[current_cohort],
           ".RData",sep=''))
print("Simulate model CHP")
m1_predict_sim <- matrix(list(),ncol=predict_sim,nrow=nrow(to_predice_obs))
for(s in 1:predict_sim){
  print(s)
  model1_par_est <- list(lambda0=sim_cohort_hp$lambda0[s],
                         eta_1=sim_cohort_hp$eta_1[s],
                         eta_2=sim_cohort_hp$eta_2[s],
                         eta_3=sim_cohort_hp$eta_3[s],
                         beta=sim_cohort_hp$beta[s],
                         f=sim_cohort_hp$f[s,])
  model1_par_matrix <- list(lambda0_matrix=matrix(model1_par_est$lambda0,
                                                  nrow=mice_number,ncol=mice_number),
                            alpha_matrix=formMatrix(function(x,y) model1_fn$alpha.fun(x,y,
                                                                                      model1_par_est$eta_1,
                                                                                      model1_par_est$eta_2,
                                                                                      model1_par_est$eta_3),
                                                    model1_par_est$f),
                            beta_matrix=matrix(model1_par_est$beta,
                                               nrow=mice_number,ncol=mice_number))
  for(cur_win in c(1:nrow(to_predice_obs))){
    m1_predict_sim[cur_win,s][[1]] <- simulateLatentHP(model1_par_matrix$lambda0_matrix,
                                                       model1_par_matrix$alpha_matrix,
                                                       model1_par_matrix$beta_matrix,
                                                       to_predice_obs$observe.time[cur_win])
  }
}

#### Predict M2 ####
load(paste(data_path,cohort_names[current_cohort],
           "/cohort_dchp_predict_stan_result_",cohort_names[current_cohort],
           ".RData",sep=''))
print("Simulate model CDCHP")
m2_predict_sim <- matrix(list(),ncol=predict_sim,nrow=nrow(to_predice_obs))
for(s in 1:predict_sim){
  print(s)
  model2_par_est <- list(gamma = sim_cohort_dchp$gamma[s,],
                         zeta = sim_cohort_dchp$zeta[s,],
                         eta_1 = sim_cohort_dchp$eta_1[s],
                         eta_2 = sim_cohort_dchp$eta_2[s],
                         eta_3 = sim_cohort_dchp$eta_3[s],
                         beta = sim_cohort_dchp$beta[s],
                         f = sim_cohort_dchp$f[s,])
  model2_par_matrix <- list(lambda0_matrix=model2_par_est$gamma%*%t(rep(1,mice_number))+
                              rep(1,mice_number)%*%t(model2_par_est$zeta),
                            alpha_matrix=formMatrix(function(x,y) model1_fn$alpha.fun(x,y,
                                                                                      model2_par_est$eta_1,
                                                                                      model2_par_est$eta_2,
                                                                                      model2_par_est$eta_3),
                                                    model2_par_est$f),
                            beta_matrix = matrix(model2_par_est$beta,
                                               nrow = mice_number,
                                               ncol = mice_number))
  for(cur_win in c(1:nrow(to_predice_obs))){
    m2_predict_sim[cur_win,s][[1]] <- simulateLatentHP(model2_par_matrix$lambda0_matrix,
                                                       model2_par_matrix$alpha_matrix,
                                                       model2_par_matrix$beta_matrix,
                                                       to_predice_obs$observe.time[cur_win])
  }
}

#### Predict M3 ####
load(paste(data_path,cohort_names[current_cohort],
           "/cohort_mmhp_predict_stan_result_",cohort_names[current_cohort],
           ".RData",sep=''))
print("Simulate model CMMHP")
m3_predict_sim <- matrix(list(),ncol = predict_sim, nrow = nrow(to_predice_obs))

for(s in 1:predict_sim){
  print(s)
  ### need to account for pair wise lambda0 and lambda1 here
  
  ### populate the matrices here
  ### need to account for pairs which aren't observed in
  ### training data
  # lambda0_matrix_curr = outer(sim_cohort_mmhp$gamma[s,],
  #                        sim_cohort_mmhp$zeta,
  #                        FUN =  "+")
  # diag(lambda0_matrix) <- 0
  
  model3_par_est <- list(lambda0 = sim_cohort_mmhp$lambda0[s],
                         lambda1 = sim_cohort_mmhp$lambda1[s],
                         eta_1=sim_cohort_mmhp$eta_1[s],
                         eta_2=sim_cohort_mmhp$eta_2[s],
                         eta_3=sim_cohort_mmhp$eta_3[s],
                         beta=sim_cohort_mmhp$beta[s],
                         f=sim_cohort_mmhp$f[s,])
  model3_par_matrix <- list(lambda0_matrix = matrix(model3_par_est$lambda0,
                                                    nrow = mice_number,
                                                    ncol = mice_number),
                              # lambda0_matrix_curr,
                            lambda1_matrix = matrix(model3_par_est$lambda1,
                                                    nrow = mice_number,
                                                    ncol = mice_number),
                              # lambda0_matrix*sim_cohort_mmhp$w_lambda[s],
                            alpha_matrix=formMatrix(function(x,y) model3_fn$alpha.fun(x,y,
                                                                                      model3_par_est$eta_1,
                                                                                      model3_par_est$eta_2),
                                                    model3_par_est$f),
                            beta_matrix=matrix(model3_par_est$beta,
                                               nrow = mice_number,
                                               ncol = mice_number),
                            q1_matrix=formMatrix(function(x,y) model3_fn$q1.fun(x,y,model3_par_est$eta_3),
                                                 model3_par_est$f),
                            q2_matrix=formMatrix(function(x,y) model3_fn$q0.fun(x,y,model3_par_est$eta_3),
                                                 model3_par_est$f))
  
  for(cur_win in c(1:nrow(to_predice_obs))){
    m3_predict_sim[cur_win,s][[1]] <- simulateLatentMMHP(lambda0_matrix = model3_par_matrix$lambda0_matrix,
                                                         lambda1_matrix = model3_par_matrix$lambda1_matrix,
                                                         alpha_matrix = model3_par_matrix$alpha_matrix,
                                                         beta_matrix = model3_par_matrix$beta_matrix,
                                                         q1_matrix = model3_par_matrix$q1_matrix,
                                                         q2_matrix = model3_par_matrix$q2_matrix,
                                                         horizon = to_predice_obs$observe.time[cur_win],
                                                         if.prefer.active = TRUE)
  }
}

#### Predict I-MMHP ####
load(paste(data_path,cohort_names[current_cohort],
           "/predict_immhp_stan_result_",cohort_names[current_cohort],
           ".RData",sep=''))
print("Simulate model IMMHP")
clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
return_df <- cleanObservationPeriod(current_cohort, clean_data)
return_df <- return_df[return_df$day<=15,]
unique_pairs_df <- return_df %>% group_by(initiator, recipient) %>%
  dplyr::summarize(count=n(),
                   observe=list(observe.id),
                   observe.length=list(observe.time),
                   no.events=list(no.events))
mmhp_predict_sim <- matrix(list(),ncol = predict_sim,
                           nrow = nrow(to_predice_obs))
mmhp_par_names <- c("lambda0","lambda1","alpha","beta","q1","q2")
mmhp_matrix_names <- c("lambda0_matrix","lambda1_matrix","alpha_matrix",
                       "beta_matrix","q1_matrix","q2_matrix")
mmhp_par_matrix <- list(lambda0_matrix = matrix(0, nrow = 12, ncol = 12),
                        lambda1_matrix = matrix(0, nrow = 12, ncol = 12),
                        alpha_matrix = matrix(0, nrow = 12, ncol = 12),
                        beta_matrix = matrix(0, nrow = 12, ncol = 12),
                        q1_matrix = matrix(0, nrow = 12, ncol = 12),
                        q2_matrix = matrix(0, nrow = 12, ncol = 12))

for(s in 1:predict_sim){
  print(s)
  for(l in c(1:length(mmhp_par_names))){
    for(pair in c(1:nrow(unique_pairs_df))){
      mmhp_par_matrix[[mmhp_matrix_names[l]]][unique_pairs_df$initiator[pair],
                                              unique_pairs_df$recipient[pair]] 
      <- sim_mmhp_sep[[mmhp_par_names[l]]][s,pair]
    }
  }
  for(cur_win in c(1:nrow(to_predice_obs))){
    mmhp_predict_sim[cur_win,s][[1]] <- simulateLatentMMHP(
      lambda0_matrix = mmhp_par_matrix$lambda0_matrix,
       lambda1_matrix = mmhp_par_matrix$lambda1_matrix,
       alpha_matrix = mmhp_par_matrix$alpha_matrix,
       beta_matrix = mmhp_par_matrix$beta_matrix,
       q1_matrix = mmhp_par_matrix$q1_matrix,
       q2_matrix = mmhp_par_matrix$q2_matrix,
       horizon = to_predice_obs$observe.time[cur_win])
  }
}

#### save this output ####

save(m1_predict_sim, m2_predict_sim, m3_predict_sim, mmhp_predict_sim,
     file = paste(data_path,cohort_names[current_cohort],
                  "/predict_simulation_",cohort_names[current_cohort],
                  ".RData",sep=''))

