# this file refits the complete C-MMHP model
# including prediction and Pearson residuals for each cohort
#

### code ###
## run this if running on the cluster
source("/rigel/stats/users/ogw2103/code/MMHP/MMHP_Latent/run_scripts/cluster_setup.R")
# ### set cohort_id based on job num
jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
cohort_id <- jobid
#cohort_id <- 1
####
save_data_path <- "output/"  #"output_june30/"

### specify the number of segments here
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
source('lib/permute_data.R')
source('lib/mmhp.R')
source('lib/uniHawkes.R')
source('lib/simulatePrediction.R')
source('lib/myGlicko.R')
source("https://gist.githubusercontent.com/jalapic/6ca3ece44bdcdc522bb735f183aa0ca0/raw/1a07f469eff08117121b6cbebcd22fb7569e3ee8/compete_extra.R")
source('lib/matrixPlotParameter.R')
source('lib/residualStructureScore.R')





#### load the data ####
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


# then fit C-MMHP here for each cohort
#######
#######
#######


#### fit the stan model ####
print(paste("Cohort",current_cohort))
stan_input_lst <- prepareData_permute_Stan(current_cohort)
### get the corresponding expert rank
permuted_data <- stan_input_lst$permuted_data

expert_rank_10_curr <- expertRankHierarchy(permuted_data)
stan_input_lst$alpha_id <- expert_rank_10_curr[1]
stan_input_lst$delta_1 <- rep(0.5,stan_input_lst$N_til)

fit_cohort_mmhp <- stan("lib/model3_current.stan",  
                        ## this will need to be updated
                        data = stan_input_lst,
                        warmup = 1000, iter = 2000, chains = 4, thin=4,
                        control=list(adapt_delta=0.999,max_treedepth=15))
sim_cohort_mmhp <- rstan::extract(fit_cohort_mmhp)
dir.create(paste(save_data_path, cohort_names[current_cohort],sep=''), recursive = TRUE, showWarnings = FALSE)
save(sim_cohort_mmhp, fit_cohort_mmhp,
     file = paste(save_data_path,cohort_names[current_cohort],
                  "/permute_cohort_mmhp_stan_result_",cohort_names[current_cohort],
                  ".RData",sep=''))