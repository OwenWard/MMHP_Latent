#### Rerun the modified simulation fits based on the modified Stan models ####
### Simulate Data from Model 3 and fit I-MMHP to each Pair

#### if running on cluster ####
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

object_par <- list(sim_lambda_1 = 0.2,
                   sim_eta_1 = 1.5, # this has to be < beta
                   gamma_var = runif(n = num_nodes, min = 0.01, max = 0.05),
                   zeta_var = runif(n = num_nodes, min = 0.01, max = 0.05),
                   sim_eta_2 = 0.6,
                   sim_eta_3 = 3,
                   sim_beta = 2,
                   f_vec_1 = seq(from = 0.05, to = 0.95,
                                 length.out = num_nodes))

object_matrix <- list(lambda0_matrix=outer(object_par$gamma_var,
                                           object_par$zeta_var, "+"),
                      lambda1_matrix=matrix(object_par$sim_lambda_1,
                                            nrow=length(object_par$f_vec_1),
                                            ncol=length(object_par$f_vec_1)),
                      alpha_matrix=formMatrix(function(x,y)
                        object_fn$alpha.fun(x,y,object_par$sim_eta_1,
                                            object_par$sim_eta_2),
                        object_par$f_vec_1),
                      beta_matrix=matrix(object_par$sim_beta,
                                         nrow=length(object_par$f_vec_1),
                                         ncol=length(object_par$f_vec_1)),
                      q1_matrix=formMatrix(function(x,y)
                        object_fn$q1.fun(x, y,
                                         object_par$sim_eta_3),
                        object_par$f_vec_1),
                      q2_matrix=formMatrix(function(x,y)
                        object_fn$q0.fun(x, y,
                                         object_par$sim_eta_3),
                        object_par$f_vec_1))


## Simulate
sim_model3_data <- list()
N_array <- array(0, c(1, num_nodes, num_nodes))
# for(i in c(1:n_sim)){
sim_model3_data <- simulateLatentMMHP(lambda0_matrix = object_matrix$lambda0_matrix,
                                      lambda1_matrix = object_matrix$lambda1_matrix,
                                      alpha_matrix = object_matrix$alpha_matrix,
                                      beta_matrix = object_matrix$beta_matrix,
                                      q1_matrix = object_matrix$q1_matrix,
                                      q2_matrix = object_matrix$q2_matrix,
                                      horizon = obs_time)
clean_sim_data <- cleanSimulationData(raw_data = sim_model3_data,
                                      cut_off = cut_off,
                                      N = length(object_par$f_vec_1))
N_array <- clean_sim_data$N_count
# }
# apply(N_array,c(2,3),mean)

dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
save(object_fn, object_par,
     object_matrix, sim_model3_data,
     file = paste(data_path,"sim_model3_immhp_", ".RData", sep=''))


##### comment out above for fitting ####
### Fit each model to simulated data ####

load(paste(data_path, "sim_model3_immhp_", ".RData", sep = ''))

clean_sim_data <- cleanSimulationData(raw_data = sim_model3_data,
                                      cut_off = cut_off,
                                      N = length(object_par$f_vec_1))


# sim_model_immhp_stan_fit <- list()
# time_vec <- clean_sim_data$event_matrix[sim_id, ]
# time_vec <- time_vec[time_vec > 0]
# diff_times <- diff(c(0, time_vec, obs_time))
# a = as.matrix(diff_times, nrow = 1, ncol = length(diff_times))
# 
# immhp_model <- cmdstan_model("lib/mmhp_single.stan")
# 
# 
# print("I-mmhp")
# ## Fit in model 3
# start_time <- Sys.time()
# sim_model_immhp_stan_fit <- immhp_model$sample(data = list(N_til = 1,
#                                                 max_Nm = length(time_vec),
#                                                 Nm = as.array(length(time_vec)),
#                                                 time_matrix= t(a),
#                                                 max_interevent =
#                                                 as.array(max(diff_times))),
#                                            chains = 4, thin = 5,
#                                            iter_sampling = 2500,
#                                            refresh = 500,
#                                            adapt_delta = 0.9)
# immhp_time <- Sys.time() - start_time
# sim_model_immhp_stan_fit$summary()
# 
# 
# sim_model_stan_sim_immhp <- sim_model_immhp_stan_fit$draws()
# post_draws <- posterior::as_draws_df(sim_model_stan_sim_immhp)
# 
# # }
# 
# ### save true parameters also
# i_ind <- clean_sim_data$I_fit[sim_id]
# j_ind <- clean_sim_data$J_fit[sim_id]
# 
# true_mmhp <- list(lambda1 = object_matrix$lambda1_matrix[i_ind, j_ind],
#                   lambda0 = object_matrix$lambda0_matrix[i_ind, j_ind],
#                   alpha = object_matrix$alpha_matrix[i_ind, j_ind],
#                   beta = object_matrix$beta_matrix[i_ind, j_ind],
#                   q1 = object_matrix$q1_matrix[i_ind, j_ind],
#                   q2 = object_matrix$q2_matrix[i_ind, j_ind])
# 
# mean_fit <- list(lambda1 = mean(post_draws$`lambda1[1]`),
#                  lambda0 = mean(post_draws$`lambda0[1]`),
#                  alpha = mean(post_draws$`alpha[1]`),
#                  beta = mean(post_draws$beta),
#                  q1 = mean(post_draws$`q1[1]`),
#                  q2 = mean(post_draws$`q2[1]`))
# 
# #### PPDiag Output
# mean_fit$Q <- matrix(c(-mean_fit$q1, mean_fit$q1, mean_fit$q2, -mean_fit$q2),
#                      nrow = 2, ncol = 2, byrow = TRUE)
# class(mean_fit) <- "mmhp"
# (resid <- pp_residual(object = mean_fit, events = time_vec, end = obs_time) )
# 
# ## Save the output
# save(true_mmhp, mean_fit,
#      resid, i_ind, j_ind,
#      sim_model_immhp_stan_fit,
#      sim_model_stan_sim_immhp,
#      file = paste(data_path,"sim_model3_fit_immhp_",
#                   sim_id,
#                   ".RData",sep=''))
# true_mmhp
# 
# cat("=========\n")
# 
# mean_fit
# 
# 
# #### Interpolate the Latent States ####
# 
# event_state_est_lst <- list()
# interpolate_state_est_lst <- list()
# 
# mmhp_par_est <- mean_fit
# 
# # delta_vec <- c(mean(sim_model_stan_sim_immhp$delta_1),
# #                1- mean(sim_model_stan_sim_immhp$delta_1))
# 
# 
# 
# viterbi_result <- myViterbi(events = time_vec,
#                             param = mmhp_par_est,
#                             termination = obs_time)
# latent_inter <- interpolateLatentTrajectory(mmhp_par_est,
#                                             time_vec,
#                                             viterbi_result$zt_v,
#                                             initial.state =
#                                               viterbi_result$initial_state,
#                                             termination.time = obs_time,
#                                             termination.state =
#                                               viterbi_result$termination_state)
# event_state_est_lst <- viterbi_result
# interpolate_state_est_lst <- latent_inter
# 
# 
# save(event_state_est_lst, interpolate_state_est_lst,
#      file = paste(data_path,
#                   "indep_mmhp_state_est_",
#                   sim_id,
#                   ".RData", sep=''))


#### Fit C-MMHP Model Instead ####

cmmhp_model <- cmdstan_model("lib/sim_model3_dc.stan")

data_list <- list(max_Nm=max(clean_sim_data$N_count),
                  N_til = length(clean_sim_data$I_fit),
                  M=sum(clean_sim_data$N_count>=cut_off),
                  N=length(object_par$f_vec_1),
                  I_fit = clean_sim_data$I_fit,
                  J_fit = clean_sim_data$J_fit,
                  T = tail(clean_sim_data$day_hour,1),
                  Nm = as.vector(clean_sim_data$N_count[clean_sim_data$N_count>=cut_off]),
                  event_matrix = clean_sim_data$event_matrix,
                  interevent_time_matrix = clean_sim_data$time_matrix,
                  max_interevent = clean_sim_data$max_interevent)


sim_cmmhp_stan_fit <- cmmhp_model$sample(data = data_list,
                                         chains = 4, thin = 5,
                                         iter_sampling = 2500,
                                         refresh = 500,
                                         adapt_delta = 0.9)

sim_model_stan_cmmhp <- sim_cmmhp_stan_fit$draws()
post_draws <- posterior::as_draws_df(sim_model_stan_cmmhp)

save(sim_model_stan_cmmhp,
     post_draws,
     file = "sim_model3_fit_cmmhp_.RData")






#### does the simulation give reasonable data compared to real data

object_par$f_vec_1
