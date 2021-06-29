#### Rerun the modified simulation fits based on the modified Stan models ####
### Simulate Data from Model 3 and fit each of the 3 models to this data

#### if running on cluster ####
source("/moto/stats/users/ogw2103/Code/MMHP_Latent/run_scripts/cluster_setup.R")

data_path <- "output/revisions/sim_m3/sim_2/"

library(cmdstanr)
library(R.utils)
library(dplyr)
library(compete)
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

object_par <- list(sim_lambda_1 = 0.6,
                   gamma_var = seq(from = 0.01, to = 0.2,
                                   length.out = num_nodes),
                   zeta_var = rep(0.1, num_nodes),
                   sim_eta_1 = 1, # from 3.5
                   sim_eta_2 = 2,#1.5, # from 2.6
                   sim_eta_3 = 3, # this seems better
                   #sim_eta_3 = 7.5,
                   sim_beta = 1.5, # from 2
                   f_vec_1 = seq(from = 0.2, to = 0.9,
                                 length.out = num_nodes))


object_matrix <- list(
  lambda0_matrix = outer(
    object_par$gamma_var,
    object_par$zeta_var, "+"
  ),
  # lambda1_matrix = matrix(object_par$sim_lambda_1,
  #   nrow = length(object_par$f_vec_1),
  #   ncol = length(object_par$f_vec_1)
  # ),
  lambda1_matrix = 1.5*lambda0_matrix,
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

dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
save(object_fn, object_par,
     object_matrix, sim_model3_data,
     file = paste(data_path,"sim_model3_",sim_id,".RData",sep=''))



### Fit each model to simulated data ####

load(paste(data_path, "sim_model3_", sim_id, ".RData", sep = ''))

## fit in the model
sim_model3_stan_fit1 <- list()
sim_model3_stan_fit2 <- list()
sim_model3_stan_fit3 <- list()

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
                  Nm = as.vector(clean_sim_data$N_count[clean_sim_data$N_count>=
                                                          cut_off]),
                  event_matrix = clean_sim_data$event_matrix,
                  interevent_time_matrix = clean_sim_data$time_matrix,
                  max_interevent = clean_sim_data$max_interevent)
print(paste(i,"model1"))
## Fit in model 1
# start_time <- Sys.time()
sim_model3_stan_fit1 <- model1$sample(data = data_list,
                                      iter_sampling = 1000,
                                      # iter_warmup = 25,
                                      refresh = 200,
                                      chains = 4)

sim_model3_fit_1 <- sim_model3_stan_fit1$draws()
sim_model3_fit1_draws <- posterior::as_draws_df(sim_model3_fit_1)
# m1_time <- Sys.time() - start_time

# print("model2")
# ## Fit in model 2
# start_time <- Sys.time()
sim_model3_stan_fit2 <- model2$sample(data = data_list,
                                      iter_sampling = 1000,
                                      refresh = 100,
                                      chains = 4)
# m2_time <- Sys.time() - start_time

sim_model3_fit_2 <- sim_model3_stan_fit2$draws()
sim_model3_fit2_draws <- posterior::as_draws_df(sim_model3_fit_2)

print("model3")
## Fit in model 3

count_data <- get_wl_matrix(df = cbind(
  clean_sim_data$start,
  clean_sim_data$end
))
isi.out <- compete::isi98(m = count_data, random = TRUE)
top_rank <- as.numeric(isi.out$best_order[1])
data_list$alpha_id <- top_rank

sim_model3_stan_fit3 <- model3$sample(data = data_list,
                                      iter_warmup = 1000,
                                      iter_sampling = 1000,
                                      chains = 4,
                                      refresh = 100)

sim_model3_fit_3 <- sim_model3_stan_fit3$draws()
sim_model3_fit3_draws <- posterior::as_draws_df(sim_model3_fit_3)

# save the summary
stan_fit3_summ <- sim_model3_stan_fit3$summary()

# }

## Extract the model fit
sim_model3_stan_sim1 <- sim_model3_fit1_draws
sim_model3_stan_sim2 <- sim_model3_fit2_draws
sim_model3_stan_sim3 <- sim_model3_fit3_draws


# m1_rank <- order(apply(sim_model3_stan_sim1$f, 2, mean))
# m2_rank <- order(apply(sim_model3_stan_sim2$f, 2, mean))
# m3_rank <- order(apply(sim_model3_stan_sim3$f, 2, mean))

## Save the output
save(object_fn, object_par, object_matrix, sim_model3_data,
     sim_model3_stan_fit1, sim_model3_stan_sim1,
     sim_model3_stan_fit2, sim_model3_stan_sim2,
     sim_model3_stan_fit3, sim_model3_stan_sim3,
     file = paste(data_path,"sim_model3_fit123_",
                  sim_id,
     ".RData",sep=''))



#### Interpolate the Latent States ####

event_state_est_lst <- list()
interpolate_state_est_lst <- list()
# for(cur_sim in c(1:n_sim)){
#   print(cur_sim)
event_state_est_lst <- matrix(list(), 
                              nrow = length(object_par$f_vec_1),
                              ncol = length(object_par$f_vec_1))
interpolate_state_est_lst <- matrix(list(),
                                    nrow = length(object_par$f_vec_1),
                                    ncol = length(object_par$f_vec_1))

lambda0_pars <- sim_model3_stan_sim3 %>% 
  select(starts_with("lambda0"))

lambda1_pars <- sim_model3_stan_sim3 %>% 
  select(starts_with("lambda1"))

f_pars <- sim_model3_stan_sim3 %>% 
  select(starts_with("f"))

###
# check the lengths match here
if(ncol(lambda0_pars) != sum(clean_sim_data$N_count >2)) {
  stop("Mismatch between params and event pairs")
}
###


lambda_0_est <- apply(lambda0_pars, 2, mean)
lambda_1_est <- apply(lambda1_pars, 2, mean)
## then put these into a matrix
lam0_par_est <- matrix(0, nrow = length(object_par$f_vec_1),
                       ncol = length(object_par$f_vec_1))
lam1_par_est <- matrix(0, nrow = length(object_par$f_vec_1),
                       ncol = length(object_par$f_vec_1))
for(i in seq_along(lambda_0_est)) {
  row_id <- clean_sim_data$I_fit[i]
  col_id <- clean_sim_data$J_fit[i]
  lam0_par_est[row_id, col_id] <- lambda_0_est[i]
  lam1_par_est[row_id, col_id] <- lambda_1_est[i]
}


mmhp_par_est <- list(lambda0 = lam0_par_est,
                     lambda1 = lam1_par_est,
                     eta_1 = mean(sim_model3_stan_sim3$eta_1),
                     eta_2 = mean(sim_model3_stan_sim3$eta_2),
                     eta_3 = mean(sim_model3_stan_sim3$eta_3),
                     beta = mean(sim_model3_stan_sim3$beta),
                     f = apply(f_pars, 2, mean))
clean_sim_data <- cleanSimulationData(raw_data = sim_model3_data, 
                                      cut_off = 1,
                                      N = length(object_par$f_vec_1))
for(cur_i in c(1:length(object_par$f_vec_1))){
  for(cur_j in c(1:length(object_par$f_vec_1))[-cur_i]){
    test.mmhp <- sim_model3_data$mmhp_matrix[cur_i, cur_j][[1]]
    object_hat <- list(lambda0 = mmhp_par_est$lambda0[cur_i, cur_j],
                       lambda1 = mmhp_par_est$lambda1[cur_i, cur_j],
                       alpha = model3_fn$alpha.fun(mmhp_par_est$f[cur_i],
                                                 mmhp_par_est$f[cur_j],
                                                 mmhp_par_est$eta_1,
                                                 mmhp_par_est$eta_2),
                       beta = mmhp_par_est$beta,
                       q1 = model3_fn$q1.fun(mmhp_par_est$f[cur_i],
                                             mmhp_par_est$f[cur_j],
                                             mmhp_par_est$eta_3),
                       q2=model3_fn$q0.fun(mmhp_par_est$f[cur_i],
                                           mmhp_par_est$f[cur_j],
                                           mmhp_par_est$eta_3))
    viterbi_result <- myViterbi(events = test.mmhp$tau[-1], 
                                param = object_hat,
                                termination = obs_time)
    latent_inter <- interpolateLatentTrajectory(object_hat, 
                                                test.mmhp$tau[-1], 
                                                viterbi_result$zt_v,
                                                initial.state = 
                                                  viterbi_result$initial_state,
                                                termination.time = obs_time,
                                                termination.state = 
                                                  viterbi_result$termination_state)
    event_state_est_lst[cur_i,cur_j][[1]] <- viterbi_result
    interpolate_state_est_lst[cur_i,cur_j][[1]] <- latent_inter
  }
}
# }

save(event_state_est_lst, interpolate_state_est_lst,
     file = paste(data_path, "fit123_state_est_", sim_id, ".RData", sep=''))


# cat("Model 1 Time:", m1_time, "\n")
# cat("Model 2 Time:", m2_time, "\n")
# cat("Model 3 Time:", m3_time, "\n")


#### Get Other Rankings

# count_data_dc <- get_wl_matrix(df = cbind(clean_sim_data$start, 
#                                           clean_sim_data$end))
# isi_dc.out <- compete::isi98(m = count_data_dc, random = TRUE)
# isi_rank_dc <- as.numeric(rev(isi_dc.out$best_order))
# 
# 
# agg_rank_data <- clean_sim_data_dc$N_count
# agg_rank_model <- stan_model(here("lib","latent_rank_agg_sim.stan"))
# 
# agg_rank_fit <- rstan::sampling(agg_rank_model,
#                                 data = list(n = 5,
#                                             n_matrix = agg_rank_data),
#                                 iter = 1000,
#                                 chains = 4)
# agg_sims <- rstan::extract(agg_rank_fit)
# 
# agg_rank <- order(apply(agg_sims$x, 2, mean))
# 
# 
# 
# ### glicko ranking also ####
# glicko_data <- tibble(start = clean_sim_data$start,
#                       end = clean_sim_data$end)
# 
# glicko_data <- glicko_data %>%
#   mutate(id = row_number(), score = 1) %>%
#   select(id, start, end, score)
# 
# gl_train <- my_glicko(glicko_data, history=TRUE, cval=2)
# 
# gl_train
# 
# gl_ranks <- order(gl_train$ratings$Rating)
# 
# ## then save these
# 
# output_rank <- tibble(truth = 1:num_nodes, 
#                       m1 = m1_rank,
#                       m2 = m2_rank,
#                       m3_dc = m3_dc_rank,
#                       isi = isi_rank_dc,
#                       agg = agg_rank,
#                       glicko = gl_ranks,
#                       sim = rep(sim_id,num_nodes))
# 
# saveRDS(output_rank, 
#         file = paste(data_path,"rank_sim",sim_id,".RDS",sep=''))
