### fit degree corrected cohort hawkes process

### code ###
# run this if running on the cluster
# source("run_scripts/cluster_setup.R")
# ### set cohort_id based on job num
# jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# jobid <- as.numeric(jobid)
# cohort_id <- jobid

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




#### then load in the data here ####

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


# then fit DC C-HP here for each cohort


#### fit the stan model ####
print(current_cohort)
stan_input_lst <- prepareDataStan(current_cohort)
stan_input_lst$scale <- 0.1
fit_cohort_dchp <- stan("lib/model2.stan",
                        data = stan_input_lst,
                        warmup = 2000, iter = 3000, thin = 4, chains = 4,
                        control=list(adapt_delta=0.95))
sim_cohort_dchp <- rstan::extract(fit_cohort_dchp)
dir.create(paste(save_data_path, cohort_names[current_cohort],sep=''), recursive = TRUE, showWarnings = FALSE)
save(sim_cohort_dchp, fit_cohort_dchp,
     file = paste(save_data_path,cohort_names[current_cohort],
                  "/cohort_dchp_stan_result_",cohort_names[current_cohort],
                  ".RData",sep=''))

### Predictions for this model ####
print(current_cohort)
stan_train_input_lst <- prepareDataStanTrain(current_cohort)
stan_train_input_lst$scale <- 0.1
fit_cohort_dchp <- stan("lib/model2.stan",
                        data = stan_train_input_lst,
                        warmup = 2000, iter = 3000, thin = 4, chains = 4,
                        control=list(adapt_delta=0.95))
sim_cohort_dchp <- rstan::extract(fit_cohort_dchp)
dir.create(paste(save_data_path, cohort_names[current_cohort],sep=''), recursive = TRUE, showWarnings = FALSE)
save(sim_cohort_dchp, fit_cohort_dchp,
     file = paste(save_data_path,cohort_names[current_cohort],
                  "/cohort_dchp_predict_stan_result_",cohort_names[current_cohort],
                  ".RData",sep=''))

#### then the pearson residuals for this fit ####
mice_number <- 12

print(current_cohort)
clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
return_df <- cleanObservationPeriod(current_cohort, clean_data)
unique_pairs_df <- return_df %>% group_by(initiator, recipient) %>%
  summarize(count=n(),
            observe=list(observe.id),
            observe.length=list(observe.time),
            no.events=list(no.events))
unique_observe_win <- unique(return_df[,c("observe.id","observe.time")])

# M2
load(paste(save_data_path,cohort_names[current_cohort],
           "/cohort_dchp_stan_result_",cohort_names[current_cohort],
           ".RData",sep=''))

model2_par_est <- list(gamma=apply(sim_cohort_dchp$gamma,2,mean),
                       zeta=apply(sim_cohort_dchp$zeta,2,mean),
                       eta_1=mean(sim_cohort_dchp$eta_1),
                       eta_2=mean(sim_cohort_dchp$eta_2),
                       eta_3=mean(sim_cohort_dchp$eta_3),
                       beta=mean(sim_cohort_dchp$beta),
                       f=apply(sim_cohort_dchp$f,2,mean))
model2_par_matrix <- list(lambda0_matrix=model2_par_est$gamma%*%t(rep(1,mice_number))+
                            rep(1,mice_number)%*%t(model2_par_est$zeta),
                          alpha_matrix=formMatrix(function(x,y) model1_fn$alpha.fun(x,y,
                                                                                    model2_par_est$eta_1,
                                                                                    model2_par_est$eta_2,
                                                                                    model2_par_est$eta_3),
                                                  model2_par_est$f),
                          beta_matrix=matrix(model2_par_est$beta,
                                             nrow=mice_number,ncol=mice_number))
m2_residual_matrix <- matrix(0,ncol=mice_number,nrow=mice_number)

for(i in 1:mice_number){
  for(j in 1:mice_number){
    pair <- which(unique_pairs_df$initiator==i&unique_pairs_df$recipient==j)
    if(length(pair)>0&(i!=j)){
      par_est <- lapply(model2_par_matrix, function(x) x[i,j])
      names(par_est) <- c("lambda0","alpha","beta")
      current_window_vec <- unique_pairs_df$observe[[pair]]
      all_residual <- 0
      for(cur in c(1:length(current_window_vec))){ 
        cur_win <- current_window_vec[cur]
        current_event_time <- return_df[return_df$initiator==i&
                                          return_df$recipient==j&
                                          return_df$observe.id==cur_win,"event.times"][[1]]
        current_obs_time <- return_df[return_df$initiator==i&
                                        return_df$recipient==j&
                                        return_df$observe.id==cur_win,"observe.time"]
        all_residual <- all_residual + uniHawkesPearsonResidual(object=par_est,
                                                                events=current_event_time,
                                                                termination = current_obs_time)
      }
      m2_residual_matrix[i,j] <- all_residual
    }
  }
}

### then save this dchp_residual_matrix ####
saveRDS(m2_residual_matrix,
        file = paste(save_data_path,cohort_names[current_cohort],
                     "/dchp_pr_matrix_",cohort_names[current_cohort],
                     ".RDS",sep=''))

