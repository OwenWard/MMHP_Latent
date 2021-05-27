##### compare inferred rank for simulated data using i&si method

source("/moto/stats/users/ogw2103/Code/MMHP_Latent/run_scripts/cluster_setup.R")


data_path <- "output/revisions/rank_sims/"


jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
sim_id <- jobid

library(here)
library(tidyverse)
library(cmdstanr)
library(compete)
library(PlayerRatings)

options(mc.cores = parallel::detectCores())

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
num_nodes <- 20

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


# #### Adjust for degree corrected sim
# 
# gamma_var <- c(0.01, 0.02, 0.03, 0.06, 0.07)
# zeta_var <- c(0.075, 0.02, 0.03, 0.05, 0.08 )
# w_lambda <- 0.4
# lambda_matrix <- outer(gamma_var, zeta_var, "+")
# # lambda_matrix
# 
# object_matrix <- object_matrix
# object_matrix$lambda0_matrix <- lambda_matrix
# 
# max_lam <- max(lambda_matrix)
# object_matrix$lambda1_matrix <- matrix(max_lam*(1 + w_lambda),
#                                           nrow = length(object_par$f_vec_1),
#                                           ncol = length(object_par$f_vec_1))


sim_model3_data <- simulateLatentMMHP(
  lambda0_matrix = object_matrix$lambda0_matrix,
  lambda1_matrix = object_matrix$lambda1_matrix,
  alpha_matrix = object_matrix$alpha_matrix,
  beta_matrix = object_matrix$beta_matrix,
  q1_matrix = object_matrix$q1_matrix,
  q2_matrix = object_matrix$q2_matrix,
  horizon = 50)
clean_sim_data <- cleanSimulationData(raw_data = sim_model3_data, 
                                         cut_off = cut_off, 
                                         N = length(object_par$f_vec_1))


### then fit the stan model (not cmdstanr)

#### model 1 ####
model1 <- cmdstan_model(here("lib", "sim_model1.stan"))

data_list <- list(max_Nm = max(clean_sim_data$N_count),
                     N_til = length(clean_sim_data$I_fit),
                     M=sum(clean_sim_data$N_count>=cut_off),
                     N=length(object_par$f_vec_1),
                     I_fit = clean_sim_data$I_fit,
                     J_fit = clean_sim_data$J_fit,
                     T = tail(clean_sim_data$day_hour,1),
                     Nm = as.vector(clean_sim_data$N_count[
                       clean_sim_data$N_count 
                       >= cut_off]),
                     event_matrix = clean_sim_data$event_matrix,
                     interevent_time_matrix = clean_sim_data$time_matrix,
                     max_interevent = clean_sim_data$max_interevent)

model1_stan_fit <- model1$sample(data = data_list,
                                 iter_sampling = 1000,
                                 chains = 4,
                                 refresh = 100)


stansims_m1 <- model1_stan_fit$summary("f")

m1_rank <- order(stansims_m1$mean)

#### model 2 ####

model2 <- cmdstan_model(here("lib", "sim_model2.stan"))

model2_stan_fit <- model2$sample(data = data_list,
                                 iter_sampling = 1000,
                                 chains = 4,
                                 refresh = 100)


stansims_m2 <- model2_stan_fit$summary("f")

m2_rank <- order(stansims_m2$mean)


#### model 3 ####
model3 <- cmdstan_model(here("lib","sim_model3.stan"))
data_list <- list(max_Nm = max(clean_sim_data$N_count),
                     N_til = length(clean_sim_data$I_fit),
                     M=sum(clean_sim_data$N_count>=cut_off),
                     N=length(object_par$f_vec_1),
                     I_fit = clean_sim_data$I_fit,
                     J_fit = clean_sim_data$J_fit,
                     T = tail(clean_sim_data$day_hour,1),
                     Nm = as.vector(clean_sim_data$N_count[
                       clean_sim_data$N_count 
                       >= cut_off]),
                     event_matrix = clean_sim_data$event_matrix,
                     interevent_time_matrix = clean_sim_data$time_matrix,
                     max_interevent = clean_sim_data$max_interevent)



model3_stan_fit <- model3$sample(data = data_list,
                                 iter_sampling = 1000,
                                 chains = 4,
                                 refresh = 100)


stansims <- model3_stan_fit$summary("f")

m3_rank <- order(stansims$mean)


### then get i&si ranking for this


count_data <- get_wl_matrix(df = cbind(clean_sim_data$start, 
                                          clean_sim_data$end))
isi.out <- compete::isi98(m = count_data, random = TRUE)
isi_rank <- as.numeric(rev(isi.out$best_order)) 


### compute agg ranking ####
agg_rank_data <- clean_sim_data$N_count
agg_rank_model <- cmdstan_model(here("lib","latent_rank_agg_sim.stan"))

agg_rank_fit <- agg_rank_model$sample(data = list(n = num_nodes,
                                                  n_matrix = agg_rank_data),
                                      iter_sampling = 1000,
                                      chains = 4)

agg_sims <- agg_rank_fit$summary("x") 

agg_rank <- order(agg_sims$mean)

### glicko ranking also ####
glicko_data <- tibble(start = clean_sim_data$start,
                      end = clean_sim_data$end)

glicko_data <- glicko_data %>%
  mutate(id = row_number(), score = 1) %>%
  select(id, start, end, score)

gl_train <- my_glicko(glicko_data, history=TRUE, cval=2)

gl_train

gl_ranks <- order(gl_train$ratings$Rating)

## then save these

output_rank <- tibble(truth = 1:20, 
                      m1 = m1_rank,
                      m2 = m2_rank,
                      m3 = m3_rank,
                      isi = isi_rank,
                      agg = agg_rank,
                      glicko = gl_ranks,
                      sim = rep(sim_id, 20))

dir.create(data_path, recursive = TRUE, showWarnings = FALSE)

saveRDS(output_rank, 
        file = paste(data_path,"rank_sim",sim_id,".RDS",sep=''))

