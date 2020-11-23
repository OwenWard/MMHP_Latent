#### Rerun the modified simulation fits based on the modified Stan models ####


#### if running on cluster ####
source("/rigel/stats/users/ogw2103/code/MMHP/MMHP_Latent/run_scripts/cluster_setup.R")

data_path <- "output/sims_m3_dc/"

library(rstan)
options(mc.cores = parallel::detectCores())
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
n_sim <- 50
cut_off <- 3

model1_fn <- list(alpha.fun = function(x,y,eta1,eta2,eta3){return(eta1*x*y*exp(-eta2*abs(x-y))/(1+exp(-eta3*(x-y))))})

model3_fn <- list(alpha.fun = function(x,y,eta1,eta2){return(eta1*x*y*exp(-eta2*abs(x-y)))},
                  q1.fun = function(x,y,eta3){return(exp(-eta3*x))},
                  q0.fun = function(x,y,eta3){return(exp(-eta3*y))})


#### Save the simulation parameters ####

object_fn <- list(alpha.fun = function(x,y,eta1,eta2){return(eta1*x*y*exp(-eta2*abs(x-y)))},
                  q1.fun = function(x,y,eta3){return(exp(-eta3*x))},
                  q0.fun = function(x,y,eta3){return(exp(-eta3*y))})

object_par <- list(sim_lambda_0 = 0.08,
                   sim_lambda_1 = 0.16,
                   sim_eta_1=2.5,
                   gamma_var = c(0.01, 0.02, 0.03, 0.06, 0.07),
                   zeta_var = c(0.075, 0.02, 0.03, 0.05, 0.08),
                   sim_eta_2 = 0.6,
                   sim_eta_3 = 5,
                   sim_beta = 1.5,
                   f_vec_1 = c(0.1, 0.2, 0.4, 0.7, 0.9))

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

# matrixPlotParameter(object_matrix$alpha_matrix)
# object_matrix$q1_matrix
# object_matrix$q2_matrix/(object_matrix$q2_matrix+object_matrix$q1_matrix)

## Simulate
# sim_model3_data <- list()
N_array <- array(0,c(1,5,5))
# for(i in c(1:n_sim)){
sim_model3_data <- simulateLatentMMHP(lambda0_matrix = object_matrix$lambda0_matrix,
                                      lambda1_matrix = object_matrix$lambda1_matrix,
                                      alpha_matrix = object_matrix$alpha_matrix,
                                      beta_matrix = object_matrix$beta_matrix,
                                      q1_matrix = object_matrix$q1_matrix,
                                      q2_matrix = object_matrix$q2_matrix,
                                      horizon = 200)
clean_sim_data <- cleanSimulationData(raw_data = sim_model3_data, 
                                      cut_off = cut_off, N = length(object_par$f_vec_1))
N_array <- clean_sim_data$N_count
# }
# apply(N_array,c(2,3),mean)

dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
save(object_fn, object_par, 
     object_matrix, sim_model3_data, 
     file = paste(data_path,"sim_model3_",sim_id,".RData",sep=''))



### Fit each model to simulatd data ####

load(paste(data_path, "sim_model3_", sim_id, ".RData", sep = ''))

## fit in the model
sim_model3_stan_fit1 <- list()
sim_model3_stan_fit2 <- list()
sim_model3_stan_fit3 <- list()

model1 <- stan_model("lib/sim_model1.stan")
model2 <- stan_model("lib/sim_model2.stan")
model3 <- stan_model("lib/sim_model3_dc.stan")

# for(i in c(1:n_sim)){
clean_sim_data <- cleanSimulationData(raw_data=sim_model3_data, 
                                      cut_off = cut_off, N = length(object_par$f_vec_1))
#stan_data <- prepareDataStan(clean_sim_data)
# these have different stan files because they're simulation models


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
print(paste(i,"model1"))
## Fit in model 1
sim_model3_stan_fit1[[1]] <- sampling(model1, data=data_list, 
                                         iter=1000, chains=4)
print("model2")
## Fit in model 2
sim_model3_stan_fit2[[1]] <- sampling(model2,
                                         data=data_list,
                                         iter=1000, chains=4)
print("model3")
## Fit in model 3
sim_model3_stan_fit3[[1]] <- sampling(model3,
                                         data=data_list,
                                         iter=1000, chains=4)
  
# }

## Extract the model fit
sim_model3_stan_sim1 <- list()
sim_model3_stan_sim2 <- list()
sim_model3_stan_sim3 <- list()
# for(i in 1:n_sim){
  sim_model3_stan_sim1 <- rstan::extract(sim_model3_stan_fit1[[1]])
  sim_model3_stan_sim2 <- rstan::extract(sim_model3_stan_fit2[[1]])
  sim_model3_stan_sim3 <- rstan::extract(sim_model3_stan_fit3[[1]])
# }

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
event_state_est_lst <- matrix(list(), nrow = length(object_par$f_vec_1),
                                      ncol = length(object_par$f_vec_1))
interpolate_state_est_lst <- matrix(list(),
                                    nrow = length(object_par$f_vec_1),
                                    ncol = length(object_par$f_vec_1))

### this doesn't take into account the degree corrected nature 
### of lambda0 and lambda1 now, which is important!

lambda_0_est <- apply(sim_model3_stan_sim3$lambda0, 2, mean)
lambda_1_est <- apply(sim_model3_stan_sim3$lambda1, 2, mean)
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
                     f = apply(sim_model3_stan_sim3$f, 2, mean))
clean_sim_data <- cleanSimulationData(raw_data=sim_model3_data, 
                                      cut_off = 1, N = length(object_par$f_vec_1))
for(cur_i in c(1:length(object_par$f_vec_1))){
  for(cur_j in c(1:length(object_par$f_vec_1))[-cur_i]){
    test.mmhp <- sim_model3_data$mmhp_matrix[cur_i,cur_j][[1]]
    object_hat <- list(lambda0=mmhp_par_est$lambda0[cur_i, cur_j],
                       lambda1=mmhp_par_est$lambda1[cur_i, cur_j],
                       alpha=model3_fn$alpha.fun(mmhp_par_est$f[cur_i],
                                                 mmhp_par_est$f[cur_j],
                                                 mmhp_par_est$eta_1,
                                                 mmhp_par_est$eta_2),
                       beta=mmhp_par_est$beta,
                       q1=model3_fn$q1.fun(mmhp_par_est$f[cur_i],
                                           mmhp_par_est$f[cur_j],
                                           mmhp_par_est$eta_3),
                       q2=model3_fn$q0.fun(mmhp_par_est$f[cur_i],
                                           mmhp_par_est$f[cur_j],
                                           mmhp_par_est$eta_3))
    viterbi_result <- myViterbi(events = test.mmhp$tau[-1], 
                                param = object_hat,
                                termination = 200)
    latent_inter <- interpolateLatentTrajectory(object_hat, 
                                                test.mmhp$tau[-1], 
                                                viterbi_result$zt_v,
                                                initial.state = viterbi_result$initial_state,
                                                termination.time = 200,
                                                termination.state = viterbi_result$termination_state)
    event_state_est_lst[[cur_sim]][cur_i,cur_j][[1]] <- viterbi_result
    interpolate_state_est_lst[[cur_sim]][cur_i,cur_j][[1]] <- latent_inter
  }
}
# }

save(event_state_est_lst, interpolate_state_est_lst,
     file = paste(data_path, "fit123_state_est_", sim_id, ".RData", sep=''))
