#### Rerun the modified simulation fits based on the modified Stan models ####

data_path <- "../output/sims/"


#### Load source code from ../lib/ ####
source('lib/uniHawkes.R')
source('lib/mmhp.R')
source('lib/simulatePrediction.R')
source('lib/plotUtil.R')
source('lib/inferLatentMmhp.R')
source('lib/drawIntensity.R')
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

object_par <- list(sim_lambda_0=0.08,
                   sim_lambda_1=0.15,
                   sim_eta_1=2.5,
                   sim_eta_2=0.6,
                   sim_eta_3=5,
                   sim_beta=1.5,
                   f_vec_1=c(0.1,0.2,0.4,0.7,0.9))

object_matrix <- list(lambda0_matrix=matrix(object_par$sim_lambda_0,
                                            nrow=length(object_par$f_vec_1),
                                            ncol=length(object_par$f_vec_1)),
                      lambda1_matrix=matrix(object_par$sim_lambda_1,
                                            nrow=length(object_par$f_vec_1),
                                            ncol=length(object_par$f_vec_1)),
                      alpha_matrix=formMatrix(function(x,y)
                        object_fn$alpha.fun(x,y,object_par$sim_eta_1,object_par$sim_eta_2),
                        object_par$f_vec_1),
                      beta_matrix=matrix(object_par$sim_beta,
                                         nrow=length(object_par$f_vec_1),
                                         ncol=length(object_par$f_vec_1)),
                      q1_matrix=formMatrix(function(x,y) object_fn$q1.fun(x,y,object_par$sim_eta_3),
                                           object_par$f_vec_1),
                      q2_matrix=formMatrix(function(x,y) object_fn$q0.fun(x,y,object_par$sim_eta_3),
                                           object_par$f_vec_1))

# matrixPlotParameter(object_matrix$alpha_matrix)
# object_matrix$q1_matrix
# object_matrix$q2_matrix/(object_matrix$q2_matrix+object_matrix$q1_matrix)

## Simulate
sim_model3_data <- rep(list(),n_sim)
N_array <- array(0,c(n_sim,5,5))
for(i in c(1:n_sim)){
  sim_model3_data[i][[1]] <- simulateLatentMMHP(lambda0_matrix=object_matrix$lambda0_matrix,
                                                lambda1_matrix=object_matrix$lambda1_matrix,
                                                alpha_matrix=object_matrix$alpha_matrix,
                                                beta_matrix=object_matrix$beta_matrix,
                                                q1_matrix=object_matrix$q1_matrix,
                                                q2_matrix=object_matrix$q2_matrix,
                                                horizon=200)
  clean_sim_data <- cleanSimulationData(raw_data=sim_model3_data[i][[1]], 
                                        cut_off = cut_off, N = length(object_par$f_vec_1))
  N_array[i,,] <- clean_sim_data$N_count
}
# apply(N_array,c(2,3),mean)

save(object_fn, object_par, 
     object_matrix, sim_model3_data, 
     file = paste(data_path,"sim_model3.RData",sep=''))