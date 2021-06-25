library(ggplot2)
library(bayesplot)
library(ggplot2)
library(posterior)
library(here)

sim_data_path <- here("output", "revisions", "sim_design", "/")

# object_par <- list(sim_lambda_1 = 0.4,
#                    gamma_var = rep(0.15, num_nodes),
#                    zeta_var = rep(0.15, num_nodes),
#                    sim_eta_1 = 3.5,
#                    sim_eta_2 = 2.6,
#                    sim_eta_3 = 7.5,
#                    sim_beta = 2,
#                    f_vec_1 = seq(from = 0.05, to = 0.95,
#                                  length.out = num_nodes))

sim_id <- 10
load(paste(sim_data_path,"sim_model3_sim123_",
                  sim_id,
                  ".RData",sep=''))

sim_fit1_draws <- sim_model3_stan_sim1
sim_fit2_draws <- sim_model3_stan_sim2
sim_fit3_draws <- sim_model3_stan_sim3


# sim_model3_stan_sim1 <- sim_fit1_draws
# sim_model3_stan_sim2 <- sim_fit2_draws
# sim_model3_stan_sim3 <- sim_fit3_draws



####
mcmc_trace(sim_fit1_draws, pars =
             c("f[1]", "f[2]", "f[3]", "f[4]", "f[5]", "f[6]", "f[7]",
               "f[8]", "f[9]", "f[10]"))
mcmc_trace(sim_fit1_draws, pars =
             c("eta_1", "eta_2", "eta_3", "beta", "lp__"))

mcmc_trace(sim_fit1_draws, pars = "eta_1") +
  geom_hline(yintercept = object_par$sim_eta_1)

mcmc_trace(sim_fit1_draws, pars = "eta_2") +
  geom_hline(yintercept = object_par$sim_eta_2)

mcmc_trace(sim_fit1_draws, pars = "eta_3") +
  geom_hline(yintercept = object_par$sim_eta_3)

mcmc_trace(sim_fit1_draws, pars = "f[1]") +
  geom_hline(yintercept = object_par$f_vec_1[1])

mcmc_trace(sim_fit1_draws, pars = "f[2]") +
  geom_hline(yintercept = object_par$f_vec_1[2])

mcmc_trace(sim_fit1_draws, pars = "f[3]") +
  geom_hline(yintercept = object_par$f_vec_1[3])

mcmc_trace(sim_fit1_draws, pars = "f[4]") +
  geom_hline(yintercept = object_par$f_vec_1[4])

mcmc_trace(sim_fit1_draws, pars = "f[5]") +
  geom_hline(yintercept = object_par$f_vec_1[5])

mcmc_trace(sim_fit1_draws, pars = "f[10]") +
  geom_hline(yintercept = object_par$f_vec_1[10])

mcmc_trace(sim_fit1_draws, pars = "beta") +
  geom_hline(yintercept = object_par$sim_beta)

####

mcmc_trace(sim_fit2_draws, pars =
             c("f[1]", "f[2]", "f[3]", "f[4]", "f[5]",
               "f[6]", "f[7]", "f[8]", "f[9]", "f[10]"))
mcmc_trace(sim_fit2_draws, pars =
             c("eta_1", "eta_2", "eta_3", "beta", "lp__"))

mcmc_trace(sim_fit2_draws, pars = "eta_1") +
  geom_hline(yintercept = object_par$sim_eta_1)

mcmc_trace(sim_fit2_draws, pars = "eta_2") +
  geom_hline(yintercept = object_par$sim_eta_2)

mcmc_trace(sim_fit2_draws, pars = "eta_3") +
  geom_hline(yintercept = object_par$sim_eta_3)

mcmc_trace(sim_fit2_draws, pars = "f[1]") +
  geom_hline(yintercept = object_par$f_vec_1[1])

mcmc_trace(sim_fit2_draws, pars = "f[2]") +
  geom_hline(yintercept = object_par$f_vec_1[2])

mcmc_trace(sim_fit2_draws, pars = "f[3]") +
  geom_hline(yintercept = object_par$f_vec_1[3])

mcmc_trace(sim_fit2_draws, pars = "f[4]") +
  geom_hline(yintercept = object_par$f_vec_1[4])

mcmc_trace(sim_fit2_draws, pars = "f[5]") +
  geom_hline(yintercept = object_par$f_vec_1[5])

mcmc_trace(sim_fit2_draws, pars = "beta") +
  geom_hline(yintercept = object_par$sim_beta)


####

mcmc_trace(sim_fit3_draws, pars =
             c("f[1]", "f[2]", "f[3]", "f[4]", "f[5]", "f[6]",
               "f[7]", "f[8]", "f[9]", "f[10]"))
mcmc_trace(sim_fit3_draws, pars =
             c("eta_1", "eta_2", "eta_3", "beta", "lp__"))

mcmc_trace(sim_fit3_draws, pars = "eta_1") +
  geom_hline(yintercept = object_par$sim_eta_1)

mcmc_trace(sim_fit3_draws, pars = "eta_2") +
  geom_hline(yintercept = object_par$sim_eta_2)

mcmc_trace(sim_fit3_draws, pars = "eta_3") +
  geom_hline(yintercept = object_par$sim_eta_3)

mcmc_trace(sim_fit3_draws, pars = "f[1]") +
  geom_hline(yintercept = object_par$f_vec_1[1])

mcmc_trace(sim_fit3_draws, pars = "f[2]") +
  geom_hline(yintercept = object_par$f_vec_1[2])

mcmc_trace(sim_fit3_draws, pars = "f[3]") +
  geom_hline(yintercept = object_par$f_vec_1[3])

mcmc_trace(sim_fit3_draws, pars = "f[4]") +
  geom_hline(yintercept = object_par$f_vec_1[4])

mcmc_trace(sim_fit3_draws, pars = "f[5]") +
  geom_hline(yintercept = object_par$f_vec_1[5])

mcmc_trace(sim_fit3_draws, pars = "f[6]") +
  geom_hline(yintercept = object_par$f_vec_1[6])

mcmc_trace(sim_fit3_draws, pars = "f[7]") +
  geom_hline(yintercept = object_par$f_vec_1[7])

mcmc_trace(sim_fit3_draws, pars = "f[8]") +
  geom_hline(yintercept = object_par$f_vec_1[8])

mcmc_trace(sim_fit3_draws, pars = "f[9]") +
  geom_hline(yintercept = object_par$f_vec_1[9])

mcmc_trace(sim_fit3_draws, pars = "f[10]") +
  geom_hline(yintercept = object_par$f_vec_1[10])

mcmc_trace(sim_fit3_draws, pars = "beta") +
  geom_hline(yintercept = object_par$sim_beta)


### compute rhats from draws, etc
rhat(extract_variable_matrix(as_draws_array(sim_fit3_draws), "eta_1"))
rhat(extract_variable_matrix(as_draws_array(sim_fit3_draws), "eta_2"))     
rhat(extract_variable_matrix(as_draws_array(sim_fit3_draws), "eta_3"))
rhat(extract_variable_matrix(as_draws_array(sim_fit3_draws), "beta"))
rhat(extract_variable_matrix(as_draws_array(sim_fit3_draws), "f[1]"))
rhat(extract_variable_matrix(as_draws_array(sim_fit3_draws), "f[2]"))
rhat(extract_variable_matrix(as_draws_array(sim_fit3_draws), "f[3]"))
rhat(extract_variable_matrix(as_draws_array(sim_fit3_draws), "f[4]"))
rhat(extract_variable_matrix(as_draws_array(sim_fit3_draws), "f[5]"))
rhat(extract_variable_matrix(as_draws_array(sim_fit3_draws), "f[6]"))
rhat(extract_variable_matrix(as_draws_array(sim_fit3_draws), "f[7]"))
rhat(extract_variable_matrix(as_draws_array(sim_fit3_draws), "f[8]"))
rhat(extract_variable_matrix(as_draws_array(sim_fit3_draws), "f[9]"))
rhat(extract_variable_matrix(as_draws_array(sim_fit3_draws), "f[10]"))

