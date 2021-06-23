### simulation
source("/moto/stats/users/ogw2103/Code/MMHP_Latent/run_scripts/cluster_setup.R")

jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

library(cmdstanr)
library(here)
library(bayesplot)
library(R.utils)
library(dplyr)
# library(compete)
options(mc.cores = parallel::detectCores())

sim_data_path <- here("output", "revisions", "sim_design", "/")


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
num_nodes <- 10
cut_off <- 3
obs_time <- 50

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

object_par <- list(sim_lambda_1 = 0.4,
                   gamma_var = rep(0.15, num_nodes),
                   zeta_var = rep(0.05, num_nodes),
                   sim_eta_1 = 3.5,
                   sim_eta_2 = 2.6,
                   sim_eta_3 = 7.5,
                   sim_beta = 2,
                   f_vec_1 = seq(from = 0.05, to = 0.95,
                                 length.out = num_nodes))

object_matrix <- list(
  lambda0_matrix = outer(
    object_par$gamma_var,
    object_par$zeta_var, "+"
  ),
  lambda1_matrix = matrix(object_par$sim_lambda_1,
                          nrow = length(object_par$f_vec_1),
                          ncol = length(object_par$f_vec_1)
  ),
  alpha_matrix = formMatrix(
    function(x, y) {
      object_fn$alpha.fun(
        x, y, object_par$sim_eta_1,
        object_par$sim_eta_2
      )
    },
    object_par$f_vec_1
  ),
  beta_matrix = matrix(object_par$sim_beta,
                       nrow = length(object_par$f_vec_1),
                       ncol = length(object_par$f_vec_1)
  ),
  q1_matrix = formMatrix(
    function(x, y) {
      object_fn$q1.fun(
        x, y,
        object_par$sim_eta_3
      )
    },
    object_par$f_vec_1
  ),
  q2_matrix = formMatrix(
    function(x, y) {
      object_fn$q0.fun(
        x, y,
        object_par$sim_eta_3
      )
    },
    object_par$f_vec_1
  )
)

# matrixPlotParameter(object_matrix$alpha_matrix)
# object_matrix$q1_matrix
# object_matrix$q2_matrix/(object_matrix$q2_matrix+object_matrix$q1_matrix)

## Simulate
# sim_model3_data <- list()
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


### Fit each model to simulated data ####


## fit in the model
model1 <- cmdstan_model("lib/sim_model1.stan")
model2 <- cmdstan_model("lib/sim_model2.stan")
model3 <- cmdstan_model("lib/sim_model3_dc.stan")

# for(i in c(1:n_sim)){
i <- sim_id
clean_sim_data <- cleanSimulationData(raw_data = sim_model3_data, 
                                      cut_off = cut_off,
                                      N = length(object_par$f_vec_1))
#stan_data <- prepareDataStan(clean_sim_data)
# these have different stan files because they're simulation models


data_list <- list(max_Nm = max(clean_sim_data$N_count),
                  N_til = length(clean_sim_data$I_fit),
                  M = sum(clean_sim_data$N_count>=cut_off),
                  N = length(object_par$f_vec_1),
                  I_fit = clean_sim_data$I_fit,
                  J_fit = clean_sim_data$J_fit,
                  T = tail(clean_sim_data$day_hour,1),
                  Nm = 
                    as.vector(clean_sim_data$N_count[clean_sim_data$N_count>=
                                                          cut_off]),
                  event_matrix = clean_sim_data$event_matrix,
                  interevent_time_matrix = clean_sim_data$time_matrix,
                  max_interevent = clean_sim_data$max_interevent)
# print(paste(i,"model1"))
## Fit in model 1
# start_time <- Sys.time()
sim_stan_fit1 <- model1$sample(data = data_list,
                                      iter_sampling = 1000,
                                      iter_warmup = 1000,
                                      refresh = 100,
                                      chains = 4)

sim_fit1 <- sim_stan_fit1$draws()
sim_fit1_draws <- posterior::as_draws_df(sim_fit1)
# # m1_time <- Sys.time() - start_time

## save the stan object
dir.create(sim_data_path, recursive = TRUE, showWarnings = FALSE)
sim_stan_fit1$save_object(file = paste0(sim_data_path,
                                        "stan_fit_1_sim_",
                                        sim_id,
                                        ".RDS"))

### look at some mcmc diagnostics
# sim_stan_fit1$summary() %>%
#   select(variable, mean, rhat) %>%
#   print(n = 30)
# mcmc_trace(sim_fit1, pars = c("eta_1", "eta_2", "eta_3"))
# mcmc_trace(sim_fit1, pars = c("f[1]", "f[2]", "f[3]", "f[4]", "f[5]"))
# mcmc_trace(sim_fit1, pars = c("lp__", "lambda0", "beta"))


# print("model2")
# ## Fit in model 2
# start_time <- Sys.time()
sim_stan_fit2 <- model2$sample(data = data_list,
                               iter_sampling = 1000,
                               iter_warmup = 1000,
                               refresh = 100,
                               chains = 4)
# m2_time <- Sys.time() - start_time

sim_fit2 <- sim_stan_fit2$draws()
sim_fit2_draws <- posterior::as_draws_df(sim_fit2)

### look at some mcmc diagnostics
# sim_model3_stan_fit2$summary() %>% 
#   select(variable, mean, rhat) %>% 
#   print(n = 42)
# 
# 
# mcmc_trace(sim_fit2, pars = c("eta_1", "eta_2", "eta_3"))
# mcmc_trace(sim_fit2, pars = c("f[1]", "f[2]", "f[3]", "f[4]", "f[5]"))
# mcmc_trace(sim_fit2, pars = c("gamma[1]",
#                               "gamma[2]",
#                               "gamma[3]",
#                               "gamma[4]",
#                               "gamma[5]"))
# mcmc_trace(sim_fit2, pars = c("zeta[1]",
#                               "zeta[2]",
#                               "zeta[3]",
#                               "zeta[4]",
#                               "zeta[5]"))
# mcmc_trace(sim_fit2, pars = c("lp__", "beta_delta", "beta"))

sim_stan_fit2$save_object(file = paste0(sim_data_path,
                                        "stan_fit_2_sim_",
                                        sim_id,
                                        ".RDS"))

print("model3")
## Fit in model 3
start_time <- Sys.time()
sim_stan_fit3 <- model3$sample(data = data_list,
                                      iter_warmup = 1000,
                                      iter_sampling = 1000,
                                      chains = 4,
                                      refresh = 100)

# m3_time <- Sys.time() - start_time
# 
sim_fit3 <- sim_stan_fit3$draws()
sim_fit3_draws <- posterior::as_draws_df(sim_fit3)

### mcmc diagnostics
# sim_stan_fit3$summary() %>%
#   select(variable, mean, rhat) %>%
#   print(n = 100)
# 
# 
# mcmc_trace(sim_fit3, pars = c("eta_1", "eta_2", "eta_3"))
# mcmc_trace(sim_fit3, pars = c("f[1]", "f[2]", "f[3]", "f[4]", "f[5]"))
# mcmc_trace(sim_fit3, pars = c("gamma[1]",
#                               "gamma[2]",
#                               "gamma[3]",
#                               "gamma[4]",
#                               "gamma[5]"))
# mcmc_trace(sim_fit3, pars = c("zeta[1]",
#                               "zeta[2]",
#                               "zeta[3]",
#                               "zeta[4]",
#                               "zeta[5]"))
# mcmc_trace(sim_fit3, pars = c("lp__", "beta_delta", "beta"))

# save the stan fit

sim_stan_fit3$save_object(file = paste0(sim_data_path,
                                        "stan_fit_3_sim_",
                                        sim_id,
                                        ".RDS"))


## Extract the model fit
sim_model3_stan_sim1 <- sim_fit1_draws
sim_model3_stan_sim2 <- sim_fit2_draws
sim_model3_stan_sim3 <- sim_fit3_draws


# m1_rank <- order(apply(sim_model3_stan_sim1$f, 2, mean))
# m2_rank <- order(apply(sim_model3_stan_sim2$f, 2, mean))
# m3_rank <- order(apply(sim_model3_stan_sim3$f, 2, mean))


## Save the output

save(object_fn, object_par, object_matrix, sim_model3_data,
     sim_model3_stan_sim1,
     sim_model3_stan_sim2,
     sim_model3_stan_sim3,
     file = paste(sim_data_path,"sim_model3_sim123_",
                  sim_id,
                  ".RData",sep=''))



####
# sim_id <- 4
# sim_stan_fit3 <- readRDS(paste0(sim_data_path,
#                                 "stan_fit_3_sim_",
#                                 sim_id,
#                                 ".RDS"))
# 
# 
# sim_stan_check$summary()
