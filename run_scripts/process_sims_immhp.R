source("/moto/stats/users/ogw2103/Code/MMHP_Latent/run_scripts/cluster_setup.R")

data_path <- "output/revisions/"

# library(rstan)
library(cmdstanr)
library(R.utils)
library(ppdiag)
# library(compete)
options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)
jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

#Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')


#### Load source code from ../lib/ ####
source('lib/uniHawkes.R')
source('lib/mmhp.R')
source('lib/simulatePrediction.R')
source('lib/plotUtil.R')
source('lib/inferLatentMmhp.R')
source('lib/drawIntensity.R')
#source('lib/prepareDataStan.R')
#source('lib/cleanData.R')
# Define global variable
n_sim <- 1
num_nodes <- 20
cut_off <- 3
obs_time <- 200

model1_fn <- list(alpha.fun = function(x, y, eta1, eta2, eta3){
  return(eta1 * x * y * exp(-eta2 * abs(x-y))/(1 + exp(-eta3 *(x-y))))})

model3_fn <- list(alpha.fun = function(x, y, eta1, eta2){
  return(eta1*x*y*exp(-eta2*abs(x-y)))},
  q1.fun = function(x,y,eta3){return(exp(-eta3*x))},
  q0.fun = function(x,y,eta3){return(exp(-eta3*y))})


#### Save the simulation parameters ####

object_fn <- list(alpha.fun = function(x,y,eta1,eta2){
  return(eta1*x*y*exp(-eta2*abs(x-y)))},
  q1.fun = function(x,y,eta3){return(exp(-eta3*x))},
  q0.fun = function(x,y,eta3){return(exp(-eta3*y))})



load(paste(data_path, "sim_model3_immhp_", ".RData", sep = ''))

clean_sim_data <- cleanSimulationData(raw_data = sim_model3_data, 
                                      cut_off = cut_off,
                                      N = length(object_par$f_vec_1))


#### then process through all the files and get the residual matrices


fit_list <- list.files(path = data_path, pattern = "sim_model3_fit_immhp_")

length(fit_list)
## should be 380
pr_matrix <- matrix(0, nrow = num_nodes, ncol = num_nodes)
rr_matrix <- matrix(0, nrow = num_nodes, ncol = num_nodes)


for(i in seq_along(fit_list)) {
  load(paste(data_path, "sim_model3_fit_immhp_",
             i,
             ".RData", sep = ''))
  pr_matrix[i_ind, j_ind] <- resid$pearson
  rr_matrix[i_ind, j_ind] <- resid$raw
}



### then plot and save these matrices
saveRDS(pr_matrix, file = here(data_path, "pr_matrix_immhp.RDS"))
saveRDS(rr_matrix, file = here(data_path, "rr_matrix_immhp.RDS"))