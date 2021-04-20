#### extract rankings and save them also


library(compete)
library(PlayerRatings)
library(rstan)
library(tidyverse)
rstan_options(auto_write = TRUE)
source(here("lib", "simulatePrediction.R"))


data_path <- "output/sims_m3_update/"


n_sim <- 50
cut_off <- 3


for(sim_id in 1:n_sim){
  load(paste(data_path,"sim_model3_",sim_id,".RData",sep=''))
  load(paste(data_path,"sim_model3_fit123_",
             sim_id,
             ".RData",sep=''))
  
  m1_rank <- order(apply(sim_model3_stan_sim1$f, 2, mean))
  m2_rank <- order(apply(sim_model3_stan_sim2$f, 2, mean))
  m3_rank <- order(apply(sim_model3_stan_sim3$f, 2, mean))
  
  ### then the rest
  clean_sim_data <- cleanSimulationData(raw_data = sim_model3_data, 
                                           cut_off = cut_off, 
                                           N = length(object_par$f_vec_1))
  
  count_data_dc <- get_wl_matrix(df = cbind(clean_sim_data$start, 
                                            clean_sim_data$end))
  isi_dc.out <- compete::isi98(m = count_data_dc, random = TRUE)
  isi_rank_dc <- as.numeric(rev(isi_dc.out$best_order))
  
  
  agg_rank_data <- clean_sim_data$N_count
  agg_rank_model <- stan_model(here("lib","latent_rank_agg_sim.stan"))
  
  agg_rank_fit <- rstan::sampling(agg_rank_model,
                                  data = list(n = 5,
                                              n_matrix = agg_rank_data),
                                  iter = 1000,
                                  chains = 4)
  agg_sims <- rstan::extract(agg_rank_fit)
  
  agg_rank <- order(apply(agg_sims$x, 2, mean))
  
  
  
  ### glicko ranking also ####
  glicko_data <- tibble(start = clean_sim_data$start,
                        end = clean_sim_data$end)
  
  glicko_data <- glicko_data %>%
    mutate(id = row_number(), score = 1) %>%
    select(id, start, end, score)
  
  gl_train <- my_glicko(glicko_data, history=TRUE, cval=2)
  
  # gl_train
  
  gl_ranks <- order(gl_train$ratings$Rating)
  
  ## then save these
  
  output_rank <- tibble(truth = 1:5, 
                        m1 = m1_rank,
                        m2 = m2_rank,
                        m3_dc = m3_rank,
                        isi = isi_rank_dc,
                        agg = agg_rank,
                        glicko = gl_ranks,
                        sim = rep(sim_id,5))
  
  saveRDS(output_rank, 
          file = paste(data_path,"rank_sim",sim_id,".RDS",sep=''))
  
}
