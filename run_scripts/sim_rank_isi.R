##### compare inferred rank for simulated data using i&si method

source("/rigel/stats/users/ogw2103/code/MMHP/MMHP_Latent/run_scripts/cluster_setup.R")


data_path <- "output/rank_sims/sims_m3_dc_isi/"


jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

library(here)
library(tidyverse)
library(rstan)
library(compete)


source(here("lib",'uniHawkes.R'))
source(here("lib","mmhp.R"))
source(here("lib",'simulatePrediction.R'))
source(here("lib",'plotUtil.R'))
source(here("lib",'inferLatentMmhp.R'))
source(here("lib",'drawIntensity.R'))
source(here("lib",'myGlicko.R'))
source(here("lib",'naiveRankHierarchy.R'))
source(here("lib",'expertRankHierarchy.R'))
cut_off <- 3

model1_fn <- list(alpha.fun = function(x, y, eta1, eta2, eta3){
  return(eta1 * x * y *exp(-eta2 * abs(x - y))/(1 + exp(-eta3 * (x - y) ) ) ) })

model3_fn <- list(alpha.fun = function(x, y, eta1, eta2){
  return(eta1 * x * y * exp(-eta2 * abs(x - y) ) )},
  q1.fun = function(x, y, eta3){
    return(exp(-eta3 * x) )},
  q0.fun = function(x, y, eta3){
    return(exp(-eta3* y) ) })


#### Save the simulation parameters ####

object_fn <- list(alpha.fun = function(x, y, eta1, eta2){
  return(eta1* x *y *exp(-eta2 * abs(x - y) ) )},
  q1.fun = function(x, y, eta3){
    return(exp(- eta3 * x) )},
  q0.fun = function(x, y, eta3){
    return(exp(- eta3 * y) )})

object_par <- list(sim_lambda_0 = 0.08,
                   sim_lambda_1 = 0.15,
                   sim_eta_1 = 2.5,
                   sim_eta_2 = 0.6,
                   sim_eta_3 = 5,
                   sim_beta = 1.5,
                   f_vec_1=c(0.1, 0.2, 0.4, 0.7, 0.9))

object_matrix <- list(lambda0_matrix = matrix(object_par$sim_lambda_0,
                                              nrow = length(object_par$f_vec_1),
                                              ncol = length(object_par$f_vec_1)),
                      lambda1_matrix = matrix(object_par$sim_lambda_1,
                                              nrow = length(object_par$f_vec_1),
                                              ncol = length(object_par$f_vec_1)),
                      alpha_matrix = formMatrix(function(x,y)
                        object_fn$alpha.fun(x, y,
                                            object_par$sim_eta_1, 
                                            object_par$sim_eta_2),
                        object_par$f_vec_1),
                      beta_matrix = matrix(object_par$sim_beta,
                                           nrow = length(object_par$f_vec_1),
                                           ncol = length(object_par$f_vec_1)),
                      q1_matrix = formMatrix(function(x, y) 
                        object_fn$q1.fun(x, y, object_par$sim_eta_3),
                        object_par$f_vec_1),
                      q2_matrix = formMatrix(function(x, y) 
                        object_fn$q0.fun(x, y, object_par$sim_eta_3),
                        object_par$f_vec_1))


#### Adjust for degree corrected sim

gamma_var <- c(0.01, 0.02, 0.03, 0.06, 0.07)
zeta_var <- c(0.075, 0.02, 0.03, 0.05, 0.08 )
lambda_dc_matrix <- outer(gamma_var, zeta_var, "+")
# lambda_dc_matrix

object_dc_matrix <- object_matrix
object_dc_matrix$lambda0_matrix <- lambda_dc_matrix



sim_model3_data_dc <- 
  simulateLatentMMHP(lambda0_matrix = object_dc_matrix$lambda0_matrix,
                     lambda1_matrix = object_dc_matrix$lambda1_matrix,
                     alpha_matrix = object_dc_matrix$alpha_matrix,
                     beta_matrix = object_dc_matrix$beta_matrix,
                     q1_matrix = object_dc_matrix$q1_matrix,
                     q2_matrix = object_dc_matrix$q2_matrix,
                     horizon = 200)
clean_sim_data_dc <- cleanSimulationData(raw_data = sim_model3_data_dc, 
                                         cut_off = cut_off, 
                                         N = length(object_par$f_vec_1))


### then fit the stan model (not cmdstanr)

model3_dc <- stan_model(here("lib","sim_model3_dc.stan"))
data_list_dc <- list(max_Nm = max(clean_sim_data_dc$N_count),
                     N_til = length(clean_sim_data_dc$I_fit),
                     M=sum(clean_sim_data_dc$N_count>=cut_off),
                     N=length(object_par$f_vec_1),
                     I_fit = clean_sim_data_dc$I_fit,
                     J_fit = clean_sim_data_dc$J_fit,
                     T = tail(clean_sim_data_dc$day_hour,1),
                     Nm = as.vector(clean_sim_data_dc$N_count[
                       clean_sim_data_dc$N_count 
                       >= cut_off]),
                     event_matrix = clean_sim_data_dc$event_matrix,
                     interevent_time_matrix = clean_sim_data_dc$time_matrix,
                     max_interevent = clean_sim_data_dc$max_interevent)



model3_dc_stan_fit <- sampling(model3_dc, data = data_list_dc,
                               iter = 1000, chains = 4)


stansims_dc <- rstan::extract(model3_dc_stan_fit)

m3_dc_rank <- order(apply(stansims_dc$f, 2, mean))


### then get i&si ranking for this


count_data_dc <- get_wl_matrix(df = cbind(clean_sim_data_dc$start, 
                                          clean_sim_data_dc$end))
isi_dc.out <- compete::isi98(m = count_data_dc, random = TRUE)
isi_rank_dc <- as.numeric(rev(isi_dc.out$best_order)) 


## then save these

output_rank <- tibble(truth = 1:5, m3_dc = m3_dc_rank,
                      isi = isi_rank_dc)

dir.create(data_path, recursive = TRUE, showWarnings = FALSE)

saveRDS(output_rank, 
        file = paste(data_path,"rank_sim",sim_id,".RDS",sep=''))
