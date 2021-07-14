### Predictions at overall level ####


### code ###
## run this if running on the cluster
source("/moto/stats/users/ogw2103/Code/MMHP_Latent/run_scripts/cluster_setup.R")
# ### set cohort_id based on job num
# jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# jobid <- as.numeric(jobid)
# cohort_id <- jobid
# cohort_id <- 1
####
data_path <- "output/revisions/"

# library(rstan)
# options(mc.cores = parallel::detectCores())
options(dplyr.summarise.inform = FALSE)

library(PlayerRatings)
library(compete)
# library(RColorBrewer)
# library(Hmisc)
library(wCorr)
# library(tidyverse)
library(dplyr)
# library(R.utils)
# library(fields)



source("lib/naiveRankHierarchy.R")
source("lib/expertRankHierarchy.R")
source("lib/cleanData.R")
source("lib/prepareDataStan.R")
source("lib/inferLatentMmhp.R")
source("lib/plotUtil.R")
source("lib/mmhp.R")
source("lib/uniHawkes.R")
source("lib/simulatePrediction.R")
source("lib/myGlicko.R")
source("https://gist.githubusercontent.com/jalapic/6ca3ece44bdcdc522bb735f183aa0ca0/raw/1a07f469eff08117121b6cbebcd22fb7569e3ee8/compete_extra.R")
source("lib/matrixPlotParameter.R")
source("lib/residualStructureScore.R")






#### then load in the data  ####

full_data <- readRDS("data/mice.RData")
# A=c9, B=c10, C=c12, D=c15, E=c16, F=c17, G=c18, H=c37, I=c38. J=c45
cohort_names <- paste("cohort", c(9, 10, 12, 15, 16, 17, 18, 37, 38, 45),
                      sep = "")
cohort_short_names <- paste("C", c(9, 10, 12, 15, 16, 17, 18, 37, 38, 45),
                            sep = "")
cut_off <- 3
mice_number <- 12


model1_fn <- list(alpha.fun = function(x, y, eta1, eta2, eta3) {
  return(eta1 * x * y * exp(-eta2 * abs(x - y)) / (1 + exp(-eta3 * (x - y))))
})

model3_fn <- list(
  alpha.fun = function(x, y, eta1, eta2) {
    return(eta1 * x * y * exp(-eta2 * abs(x - y)))
  },
  q1.fun = function(x, y, eta3) {
    return(exp(-eta3 * x))
  },
  q0.fun = function(x, y, eta3) {
    return(exp(-eta3 * y))
  }
)

# Define the cohorts will be fitted
fit_cohorts <- c(1:10)
naive_rank_10 <- list()
expert_rank_10 <- list()
for (current_cohort in fit_cohorts) {
  naive_rank_10[[current_cohort]] <-
    naiveRankHierarchy(full_data[[cohort_names[current_cohort]]])
  expert_rank_10[[current_cohort]] <-
    expertRankHierarchy(full_data[[cohort_names[current_cohort]]])
}

predict_sim <- 1000

f_result_list <- list()
f_sd_list <- list()
for (current_cohort in 1:10) {
  f_result_list[[current_cohort]] <- matrix(0, nrow = 4, ncol = mice_number)
  f_sd_list[[current_cohort]] <- matrix(0, nrow = 4, ncol = mice_number)
  ## agg
  load(paste(data_path, cohort_names[current_cohort],
    "/agg_rank_stan_result_", cohort_names[current_cohort],
    ".RData",
    sep = ""
  ))
  x_draws <- sim_agg_rank %>% select(starts_with("x"))

  f_result_list[[current_cohort]][1, ] <- apply(
    x_draws,
    2, median
  )
  f_sd_list[[current_cohort]][1, ] <- apply(
    x_draws,
    2, sd
  )
  ## M1
  load(paste(data_path, cohort_names[current_cohort],
    "/cohort_hp_stan_result_", cohort_names[current_cohort],
    ".RData",
    sep = ""
  ))
  f_draws <- sim_cohort_hp %>% select(starts_with("f"))
  f_result_list[[current_cohort]][2, ] <- apply(f_draws, 2, median)
  f_sd_list[[current_cohort]][2, ] <- apply(f_draws, 2, sd)
  ## M2
  load(paste(data_path,
    cohort_names[current_cohort],
    "/cohort_dchp_stan_result_",
    cohort_names[current_cohort],
    ".RData",
    sep = ""
  ))
  f_draws <- sim_cohort_dchp %>% select(starts_with("f"))
  f_result_list[[current_cohort]][3, ] <- apply(f_draws, 2, median)
  f_sd_list[[current_cohort]][3, ] <- apply(f_draws, 2, sd)
  ## M3
  load(paste(data_path,
    cohort_names[current_cohort],
    "/cohort_mmhp_stan_result_",
    cohort_names[current_cohort],
    ".RData",
    sep = ""
  ))
  f_draws <- sim_cohort_mmhp %>% select(starts_with("f"))
  f_result_list[[current_cohort]][4, ] <- apply(f_draws, 2, median)
  f_sd_list[[current_cohort]][4, ] <- apply(f_draws, 2, sd)
}

spearman_df <- matrix(0, nrow = 10, ncol = 4)
pearson_df <- matrix(0, nrow = 10, ncol = 4)
for (current_cohort in fit_cohorts) {
  clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
  expert_pair <- mice_number + 1 - order(expert_rank_10[[current_cohort]])
  for (current_method in c(1:4)) {
    spearman_df[current_cohort, current_method] <- weightedCorr(expert_pair,
      f_result_list[[current_cohort]][current_method, ],
      method = "Spearman",
      weights = 1 / f_sd_list[[current_cohort]][current_method, ]
    )
    pearson_df[current_cohort, current_method] <- weightedCorr(expert_pair,
      f_result_list[[current_cohort]][current_method, ],
      method = "Pearson",
      weights = 1 / f_sd_list[[current_cohort]][current_method, ]
    )
  }
}

save(f_result_list,
  f_sd_list,
  spearman_df,
  pearson_df,
  file = paste(data_path, "weighted_rank_data.RData", sep = "")
)


#### then n matrix for all cohorts

no_method <- 5
all_cohort_norm_df <- data.frame(
  norm = rep(NA, 10 * 6 * no_method),
  day = rep(NA, 10 * 6 * no_method),
  method = rep(NA, 10 * 6 * no_method),
  cohort = rep(NA, 10 * 6 * no_method)
)

all_cohort_mae_df <- data.frame(
  norm = rep(NA, 10 * 6 * no_method),
  day = rep(NA, 10 * 6 * no_method),
  method = rep(NA, 10 * 6 * no_method),
  cohort = rep(NA, 10 * 6 * no_method)
)

predict_day_norm_df_lst <- list()
predict_day_mae_df_lst <- list()
cur_all <- 1
for (current_cohort in 1:10) {
  print(current_cohort)
  load(paste(data_path, cohort_names[current_cohort],
    "/predict_simulation_", cohort_names[current_cohort],
    ".RData",
    sep = ""
  ))
  load(paste(data_path, cohort_names[current_cohort],
    "/predict_dsnl_poisson_stan_result_", cohort_names[current_cohort],
    ".RData",
    sep = ""
  ))
  clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
  train_clean_data <- cleanDataForTraining(
    all_raw_data = full_data[[cohort_names[current_cohort]]],
    cut_off = 0, N = 12, train_period = 15
  )
  test_day_end <- unlist(lapply(c(16:21), function(x) clean_data$day_hour[max(which(clean_data$day == x))]))
  test_start <- tail(train_clean_data$day_hour, 1)
  ### this causing issue
  return_df <- cleanObservationPeriod(current_cohort,
    raw_df = full_data[[cohort_names[current_cohort]]],
    clean_data
  )
  to_predice_obs <- unique(return_df[
    return_df$day >= 16,
    c("day", "observe.id", "observe.time")
  ])

  real_N_matrix_list <- list()
  predict_day_norm_df <- data.frame(
    norm = rep(NA, predict_sim * 6 * no_method),
    day = rep(NA, predict_sim * 6 * no_method),
    method = rep(NA, predict_sim * 6 * no_method)
  )
  predict_day_mae_df <- data.frame(
    mae = rep(NA, predict_sim * 6 * no_method),
    day = rep(NA, predict_sim * 6 * no_method),
    method = rep(NA, predict_sim * 6 * no_method)
  )
  est_array_m1 <- array(0, dim = c(6, predict_sim, mice_number, mice_number))
  est_array_m2 <- array(0, dim = c(6, predict_sim, mice_number, mice_number))
  est_array_m3 <- array(0, dim = c(6, predict_sim, mice_number, mice_number))
  est_array_mmhp <- array(0, dim = c(6, predict_sim, mice_number, mice_number))
  est_array_dsnl <- array(0, dim = c(6, predict_sim, mice_number, mice_number))
  cur <- 1
  for (d_test in c(16:21)) {
    indicate_day <- (return_df$day > 15 & return_df$day <= d_test)
    temp_real <- cleanSimulationDataForNCount(list(
      start =
        return_df$initiator[indicate_day],
      end =
        return_df$recipient[indicate_day],
      day_hour =
        c(1:sum(indicate_day))
    ))$N_count
    real_N_matrix_list[[d_test]] <- temp_real

    observe_windows <- which(to_predice_obs$day == d_test)
    for (s in c(1:predict_sim)) {
      # M1
      if (d_test > 16) {
        est_array_m1[d_test - 15, s, , ] <- est_array_m1[d_test - 16, s, , ]
      }
      for (win in observe_windows) {
        est_array_m1[d_test - 15, s, , ] <- est_array_m1[d_test - 15, s, , ] +
          cleanSimulationDataForNCount(m1_predict_sim[win, s][[1]])$N_count
      }
      predict_day_norm_df[cur, "norm"] <- sqrt(sum((real_N_matrix_list[[d_test]] -
        est_array_m1[d_test - 15, s, , ])^2))
      ### try MAE for this
      predict_day_mae_df[cur, "mae"] <- mean(abs(real_N_matrix_list[[d_test]] -
        est_array_m1[d_test - 15, s, , ]))

      # M2
      if (d_test > 16) {
        est_array_m2[d_test - 15, s, , ] <- est_array_m2[d_test - 16, s, , ]
      }
      for (win in observe_windows) {
        est_array_m2[d_test - 15, s, , ] <- est_array_m2[d_test - 15, s, , ] +
          cleanSimulationDataForNCount(m2_predict_sim[win, s][[1]])$N_count
      }
      predict_day_norm_df[cur + 1, "norm"] <- sqrt(sum((real_N_matrix_list[[d_test]] -
        est_array_m2[d_test - 15, s, , ])^2))
      ### try MAE
      predict_day_mae_df[cur + 1, "mae"] <- mean(abs(real_N_matrix_list[[d_test]] -
        est_array_m2[d_test - 15, s, , ]))

      # M3
      if (d_test > 16) {
        est_array_m3[d_test - 15, s, , ] <- est_array_m3[d_test - 16, s, , ]
      }
      for (win in observe_windows) {
        est_array_m3[d_test - 15, s, , ] <- est_array_m3[d_test - 15, s, , ] +
          cleanSimulationDataForNCount(m3_predict_sim[win, s][[1]])$N_count
      }
      predict_day_norm_df[cur + 2, "norm"] <- sqrt(sum((real_N_matrix_list[[d_test]] -
        est_array_m3[d_test - 15, s, , ])^2))
      ### MAE here
      predict_day_mae_df[cur + 2, "mae"] <- mean(abs(real_N_matrix_list[[d_test]] -
        est_array_m3[d_test - 15, s, , ]))

      # I-MMHP
      if (d_test > 16) {
        est_array_mmhp[d_test - 15, s, , ] <- est_array_mmhp[d_test - 16, s, , ]
      }
      for (win in observe_windows) {
        est_array_mmhp[d_test - 15, s, , ] <- est_array_mmhp[d_test - 15, s, , ] +
          cleanSimulationDataForNCount(mmhp_predict_sim[win, s][[1]])$N_count
      }
      predict_day_norm_df[cur + 3, "norm"] <- sqrt(sum((real_N_matrix_list[[d_test]] -
        est_array_mmhp[d_test - 15, s, , ])^2))

      ### MAE
      predict_day_mae_df[cur + 3, "mae"] <- mean(abs(real_N_matrix_list[[d_test]] -
        est_array_mmhp[d_test - 15, s, , ]))

      # DSNL
      dsnl_lambda <- sim_predict_dsnl %>% select(starts_with("lambda_d"))
      dsnl_lambda_arr <- array(as.matrix(dsnl_lambda), dim = c(1000, 6, 12, 12))
      if (d_test == 16) {
        ## d_test - 15 corresponds to the first entry
        est_array_dsnl[d_test - 15, s, , ] <- dsnl_lambda_arr[s, d_test - 15, , ] *
          (test_day_end[d_test - 15] - test_start)
      } else {
        est_array_dsnl[d_test - 15, s, , ] <- est_array_dsnl[d_test - 16, s, , ] +
          dsnl_lambda_arr[s, d_test - 15, , ] *
            (test_day_end[d_test - 15] - test_day_end[d_test - 16])
      }
      predict_day_norm_df[cur + 4, "norm"] <- sqrt(sum((real_N_matrix_list[[d_test]] -
        est_array_dsnl[d_test - 15, s, , ])^2))

      ### MAE
      predict_day_mae_df[cur + 4, "mae"] <- mean(abs(real_N_matrix_list[[d_test]] -
        est_array_dsnl[d_test - 15, s, , ]))

      predict_day_norm_df[c(cur:(cur + no_method - 1)), "day"] <- rep(d_test,
                                                                      no_method)
      predict_day_norm_df[c(cur:(cur + no_method - 1)), "method"] <- c("m1",
                                                                       "m2",
                                                                       "m3",
                                                                       "mmhp",
                                                                       "dsnl")
      ### MAE
      predict_day_mae_df[c(cur:(cur + no_method - 1)), "day"] <- rep(d_test,
                                                                     no_method)
      predict_day_mae_df[c(cur:(cur + no_method - 1)), "method"] <- c("m1",
                                                                      "m2",
                                                                      "m3",
                                                                      "mmhp",
                                                                      "dsnl")
      cur <- cur + no_method
    }
    all_cohort_norm_df[cur_all, "norm"] <- sqrt(sum((real_N_matrix_list[[d_test]] -
      apply(est_array_m1[d_test - 15, , , ], c(2, 3), median))^2))
    all_cohort_norm_df[cur_all + 1, "norm"] <- sqrt(sum((real_N_matrix_list[[d_test]] -
      apply(est_array_m2[d_test - 15, , , ], c(2, 3), median))^2))
    all_cohort_norm_df[cur_all + 2, "norm"] <- sqrt(sum((real_N_matrix_list[[d_test]] -
      apply(est_array_m3[d_test - 15, , , ], c(2, 3), median))^2))
    all_cohort_norm_df[cur_all + 3, "norm"] <- sqrt(sum((real_N_matrix_list[[d_test]] -
      apply(est_array_mmhp[d_test - 15, , , ], c(2, 3), median))^2))
    all_cohort_norm_df[cur_all + 4, "norm"] <- sqrt(sum((real_N_matrix_list[[d_test]] -
      apply(est_array_dsnl[d_test - 15, , , ], c(2, 3), median))^2))
    all_cohort_norm_df[c(cur_all:(cur_all + no_method - 1)),
                       "day"] <- rep(d_test, no_method)
    all_cohort_norm_df[c(cur_all:(cur_all + no_method - 1)),
                       "cohort"] <- rep(current_cohort, no_method)
    all_cohort_norm_df[c(cur_all:(cur_all + no_method - 1)),
                       "method"] <- c("m1", "m2", "m3", "mmhp", "dsnl")
    ### repeat for mae
    all_cohort_mae_df[cur_all, "norm"] <- mean(abs(real_N_matrix_list[[d_test]] -
      apply(est_array_m1[d_test - 15, , , ], c(2, 3), median)))
    all_cohort_mae_df[cur_all + 1, "norm"] <- mean(abs(real_N_matrix_list[[d_test]] -
      apply(est_array_m2[d_test - 15, , , ], c(2, 3), median)))
    all_cohort_mae_df[cur_all + 2, "norm"] <- mean(abs(real_N_matrix_list[[d_test]] -
      apply(est_array_m3[d_test - 15, , , ], c(2, 3), median)))
    all_cohort_mae_df[cur_all + 3, "norm"] <- mean(abs(real_N_matrix_list[[d_test]] -
      apply(est_array_mmhp[d_test - 15, , , ], c(2, 3), median)))
    all_cohort_mae_df[cur_all + 4, "norm"] <- mean(abs(real_N_matrix_list[[d_test]] -
      apply(est_array_dsnl[d_test - 15, , , ], c(2, 3), median)))
    all_cohort_mae_df[c(cur_all:(cur_all + no_method - 1)), "day"] <- rep(d_test, no_method)
    all_cohort_mae_df[c(cur_all:(cur_all + no_method - 1)), "cohort"] <- rep(current_cohort, no_method)
    all_cohort_mae_df[c(cur_all:(cur_all + no_method - 1)), "method"] <- c("m1", "m2", "m3", "mmhp", "dsnl")


    cur_all <- cur_all + no_method
  }
  predict_day_norm_df_lst[[current_cohort]] <- predict_day_norm_df
  predict_day_mae_df_lst[[current_cohort]] <- predict_day_mae_df
}
save(predict_day_norm_df_lst,
  predict_day_mae_df_lst,
  all_cohort_norm_df,
  all_cohort_mae_df,
  est_array_m1,
  est_array_m2,
  est_array_m3,
  est_array_mmhp,
  est_array_dsnl,
  file = paste(data_path, "plot_N_predict_revision.RData", sep = "")
)




#### Rank Prediction

predict_sim <- 1000
no_method <- 5

all_cohort_rank_df <- data.frame(
  spearman = rep(NA, 10 * 6 * no_method),
  day = rep(NA, 10 * 6 * no_method),
  method = rep(NA, 10 * 6 * no_method),
  cohort = rep(NA, 10 * 6 * no_method)
)
predict_day_rank_df_lst <- list()
cur_all <- 1
for (current_cohort in fit_cohorts) {
  print(current_cohort)
  load(paste(data_path, cohort_names[current_cohort],
    "/predict_simulation_", cohort_names[current_cohort],
    ".RData",
    sep = ""
  ))
  load(paste(data_path, cohort_names[current_cohort],
    "/predict_dsnl_poisson_stan_result_", cohort_names[current_cohort],
    ".RData",
    sep = ""
  ))
  train_clean_data <- cleanDataForTraining(
    all_raw_data = full_data[[cohort_names[current_cohort]]],
    cut_off = 0, N = 12, train_period = 15
  )
  clean_data_all <- cleanData(full_data[[cohort_names[current_cohort]]],
    cut_off = 1
  )
  clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
  return_df <- cleanObservationPeriod(current_cohort,
    raw_df = full_data[[cohort_names[current_cohort]]],
    clean_data
  )
  to_predice_obs <- unique(return_df[(return_df$day >= 16) &
                                       (return_df$day <= 21), c("day",
                                                                "observe.id",
                                                                "observe.time")]) ## change with day info

  #------------ gl score
  df <- full_data[[cohort_names[current_cohort]]]
  df <- df[df$Actor != "Start" & df$Actor != "End", ]
  df1 <- expandrows(df)
  df1 <- df1[order(df1$Timestamp), ] # ensure in date order
  df1$event <- 1:nrow(df1)
  glick.df <- df1[,
                  c("event", "Actor", "Recipient", "score"),
                  with = FALSE] # need event, actor, recipient, score
  real.gl <- my_glicko(glick.df, history = TRUE, cval = 2)

  test_day_end <- unlist(lapply(c(16:21),
                                function(x) clean_data_all$day_hour[max(which(clean_data_all$day == x))]))
  test_day_end_idx <- unlist(lapply(c(16:21),
                                    function(x) max(which(clean_data_all$day == x))))
  test_start <- tail(train_clean_data$day_hour, 1)

  predict_day_rank_df <- data.frame(
    spearman = rep(NA, predict_sim * 6 * no_method),
    day = rep(NA, predict_sim * 6 * no_method),
    method = rep(NA, predict_sim * 6 * no_method)
  )
  est_intensity_array_m1 <- array(0,
                                  dim = c(6,
                                          predict_sim,
                                          mice_number,
                                          mice_number))
  est_intensity_array_m2 <- array(0,
                                  dim = c(6,
                                          predict_sim,
                                          mice_number,
                                          mice_number))
  est_intensity_array_m3 <- array(0,
                                  dim = c(6,
                                          predict_sim,
                                          mice_number,
                                          mice_number))
  est_intensity_array_mmhp <- array(0,
                                    dim = c(6,
                                            predict_sim,
                                            mice_number,
                                            mice_number))
  est_intensity_array_dsnl <- array(0,
                                    dim = c(6,
                                            predict_sim,
                                            mice_number,
                                            mice_number))

  # Load parameter values
  ## m1
  load(paste(data_path, cohort_names[current_cohort],
    "/cohort_hp_predict_stan_result_", cohort_names[current_cohort],
    ".RData",
    sep = ""
  ))
  ## m2
  load(paste(data_path,
    cohort_names[current_cohort],
    "/cohort_dchp_predict_stan_result_",
    cohort_names[current_cohort],
    ".RData",
    sep = ""
  ))
  ## m3
  load(paste(data_path,
    cohort_names[current_cohort],
    "/cohort_mmhp_predict_stan_result_",
    cohort_names[current_cohort],
    ".RData",
    sep = ""
  ))
  ## I-mmhp
  load(paste(data_path,
    cohort_names[current_cohort],
    "/predict_immhp_stan_result_",
    cohort_names[current_cohort],
    ".RData",
    sep = ""
  ))
  cur <- 1
  for (d_test in c(16:21)) {
    print(d_test)
    if (test_day_end_idx[d_test - 15] > dim(real.gl$history)[2]) {
      real_gl_vec <- real.gl$history[, dim(real.gl$history)[2], 1]
    } else {
      real_gl_vec <- real.gl$history[, test_day_end_idx[d_test - 15], 1]
    }

    for (s in c(1:predict_sim)) {
      ## m1
      f_draws <- sim_cohort_hp %>% select(starts_with("f"))
      model1_par_est <- list(
        lambda0 = sim_cohort_hp$lambda0[s],
        eta_1 = sim_cohort_hp$eta_1[s],
        eta_2 = sim_cohort_hp$eta_2[s],
        eta_3 = sim_cohort_hp$eta_3[s],
        beta = sim_cohort_hp$beta[s],
        f = as.numeric(f_draws[s, ])
      )
      model1_par_matrix <- list(
        lambda0_matrix = matrix(model1_par_est$lambda0,
          nrow = mice_number, ncol = mice_number
        ),
        alpha_matrix = formMatrix(
          function(x, y) {
            model1_fn$alpha.fun(
              x, y,
              model1_par_est$eta_1,
              model1_par_est$eta_2,
              model1_par_est$eta_3
            )
          },
          model1_par_est$f
        ),
        beta_matrix = matrix(model1_par_est$beta,
          nrow = mice_number, ncol = mice_number
        )
      )
      ## m2
      f_draws <- sim_cohort_dchp %>% select(starts_with("f"))
      gamma_draws <- sim_cohort_dchp %>% select(starts_with("gamma"))
      zeta_draws <- sim_cohort_dchp %>% select(starts_with("zeta"))
      model2_par_est <- list(
        gamma = as.numeric(gamma_draws[s, ]),
        zeta = as.numeric(zeta_draws[s, ]),
        eta_1 = sim_cohort_dchp$eta_1[s],
        eta_2 = sim_cohort_dchp$eta_2[s],
        eta_3 = sim_cohort_dchp$eta_3[s],
        beta = sim_cohort_dchp$beta[s],
        f = as.numeric(f_draws[s, ])
      )
      model2_par_matrix <- list(
        lambda0_matrix =
          model2_par_est$gamma %*% t(rep(1, mice_number)) +
            rep(1, mice_number) %*% t(model2_par_est$zeta),
        alpha_matrix = formMatrix(
          function(x, y) {
            model1_fn$alpha.fun(
              x, y,
              model2_par_est$eta_1,
              model2_par_est$eta_2,
              model2_par_est$eta_3
            )
          },
          model2_par_est$f
        ),
        beta_matrix = matrix(model2_par_est$beta,
          nrow = mice_number,
          ncol = mice_number
        )
      )

      ## m3
      ### need to update the matrix fill in here
      lam0_draws <- sim_cohort_mmhp %>% select(starts_with("lambda0"))
      lam1_draws <- sim_cohort_mmhp %>% select(starts_with("lambda1"))
      f_draws <- sim_cohort_mmhp %>% select(starts_with("f"))
      model3_par_est <- list(
        lambda0 = as.numeric(lam0_draws[s, ]),
        lambda1 = as.numeric(lam1_draws[s, ]),
        eta_1 = sim_cohort_mmhp$eta_1[s],
        eta_2 = sim_cohort_mmhp$eta_2[s],
        eta_3 = sim_cohort_mmhp$eta_3[s],
        beta = sim_cohort_mmhp$beta[s],
        f = as.numeric(f_draws[s, ])
      )
      ## update matrix here
      lam0_matrix <- matrix(0, nrow = mice_number, ncol = mice_number)
      lam1_matrix <- matrix(0, nrow = mice_number, ncol = mice_number)
      for (i in seq_along(clean_data$I_fit)) {
        row_id <- clean_data$I_fit[i]
        col_id <- clean_data$J_fit[i]
        lam0_matrix[row_id, col_id] <- model3_par_est$lambda0[i]
        lam1_matrix[row_id, col_id] <- model3_par_est$lambda1[i]
      }

      model3_par_matrix <- list(
        lambda0_matrix = lam0_matrix,
        lambda1_matrix = lam1_matrix,
        alpha_matrix = formMatrix(
          function(x, y) {
            model3_fn$alpha.fun(
              x, y,
              model3_par_est$eta_1,
              model3_par_est$eta_2
            )
          },
          model3_par_est$f
        ),
        beta_matrix = matrix(model3_par_est$beta,
          nrow = mice_number, ncol = mice_number
        ),
        q1_matrix = formMatrix(
          function(x, y) model3_fn$q1.fun(x, y, model3_par_est$eta_3),
          model3_par_est$f
        ),
        q2_matrix = formMatrix(
          function(x, y) model3_fn$q0.fun(x, y, model3_par_est$eta_3),
          model3_par_est$f
        )
      )

      ## mmhp
      train_return_df <- return_df[return_df$day <= 15, ]
      unique_pairs_df <- train_return_df %>%
        group_by(initiator, recipient) %>%
        dplyr::summarize(
          count = n(),
          observe = list(observe.id),
          observe.length = list(observe.time),
          no.events = list(no.events)
        )
      mmhp_par_names <- c("lambda0", "lambda1", "alpha", "beta", "q1", "q2")
      mmhp_matrix_names <- c(
        "lambda0_matrix", "lambda1_matrix", "alpha_matrix",
        "beta_matrix", "q1_matrix", "q2_matrix"
      )
      mmhp_par_matrix <- list(
        lambda0_matrix = matrix(0, nrow = 12, ncol = 12),
        lambda1_matrix = matrix(0, nrow = 12, ncol = 12),
        alpha_matrix = matrix(0, nrow = 12, ncol = 12),
        beta_matrix = matrix(0, nrow = 12, ncol = 12),
        q1_matrix = matrix(0, nrow = 12, ncol = 12),
        q2_matrix = matrix(0, nrow = 12, ncol = 12)
      )
      
      for (l in c(1:length(mmhp_par_names))) {
        for (pair in c(1:nrow(unique_pairs_df))) {
          curr_draws <- sim_mmhp_sep %>% 
            select(starts_with(mmhp_par_names[l])) %>% 
            select(!contains("delta"))
          
          mmhp_par_matrix[[mmhp_matrix_names[l]]][
            unique_pairs_df$initiator[pair],
            unique_pairs_df$recipient[pair]
          ] <- as.numeric(curr_draws[s, pair])
          # mmhp_par_matrix[[mmhp_matrix_names[l]]][
          #   unique_pairs_df$initiator[pair],
          #   unique_pairs_df$recipient[pair]
          # ] <- sim_mmhp_sep[[mmhp_par_names[l]]][s, pair]
        }
      }

      cur_day_windows <- which(to_predice_obs$day == d_test)
      for (w in cur_day_windows) { 
        # the window is absolute not change with the day
        ## M1
        clean_sim <- cleanSimulationDataForNCount(m1_predict_sim[w, s][[1]])
        for (i in c(1:mice_number)) {
          for (j in c(1:mice_number)[-i]) {
            est_par <- lapply(model1_par_matrix, function(x) x[i, j])
            names(est_par) <- c("lambda0", "alpha", "beta")
            est_intensity_array_m1[d_test - 15, s, i, j] <- max(
              est_intensity_array_m1[d_test - 15, s, i, j],
              uniHawkesIntensity(
                object = est_par,
                events = unlist(clean_sim$time_each_pair[i, j]),
                current_time = to_predice_obs$observe.time[w]
              )
            )
          }
        }

        # M2
        clean_sim <- cleanSimulationDataForNCount(m2_predict_sim[w, s][[1]])
        for (i in c(1:mice_number)) {
          for (j in c(1:mice_number)[-i]) {
            est_par <- lapply(model2_par_matrix, function(x) x[i, j])
            names(est_par) <- c("lambda0", "alpha", "beta")
            est_intensity_array_m2[d_test - 15, s, i, j] <- max(
              est_intensity_array_m2[d_test - 15, s, i, j],
              uniHawkesIntensity(
                object = est_par,
                events = unlist(clean_sim$time_each_pair[i, j]),
                current_time = to_predice_obs$observe.time[w]
              )
            )
          }
        }

        # M3
        clean_sim <- cleanSimulationDataForNCount(mmhp_predict_sim[w, s][[1]])
        for (i in c(1:mice_number)) {
          for (j in c(1:mice_number)[-i]) {
            est_par <- lapply(mmhp_par_matrix, function(x) x[i, j])
            names(est_par) <- c("lambda0",
                                "lambda1",
                                "alpha",
                                "beta",
                                "q1",
                                "q2")
            est_intensity_array_m3[d_test - 15, s, i, j] <- max(
              est_intensity_array_m3[d_test - 15, s, i, j],
              mmhpIntensityAtTime(
                params = est_par,
                events = unlist(clean_sim$time_each_pair[i, j]),
                current_time = to_predice_obs$observe.time[w],
                latent_z = 1,
                latent_x = 0
              )
            )
          }
        }

        # MMHP
        clean_sim <- cleanSimulationDataForNCount(m3_predict_sim[w, s][[1]])
        for (i in c(1:mice_number)) {
          for (j in c(1:mice_number)[-i]) {
            est_par <- lapply(model3_par_matrix, function(x) x[i, j])
            names(est_par) <- c("lambda0",
                                "lambda1",
                                "alpha",
                                "beta",
                                "q1",
                                "q2")
            est_intensity_array_mmhp[d_test - 15, s, i, j] <- max(
              est_intensity_array_mmhp[d_test - 15, s, i, j],
              mmhpIntensityAtTime(
                params = est_par,
                events = unlist(clean_sim$time_each_pair[i, j]),
                current_time = to_predice_obs$observe.time[w],
                latent_z = mmhp_predict_sim[w, s][[1]][[4]][i, j][[1]]$z,
                latent_x = mmhp_predict_sim[w, s][[1]][[4]][i, j][[1]]$x
              )
            )
          }
        }
      }

      # m1
      predict_day_rank_df[cur, "spearman"] <- cor.test(real_gl_vec,
        rowSums(est_intensity_array_m1[d_test - 15, s, , ]),
        method = "spearman", exact = FALSE
      )$estimate
      # m2
      predict_day_rank_df[cur + 1, "spearman"] <- cor.test(real_gl_vec,
        rowSums(est_intensity_array_m2[d_test - 15, s, , ]),
        method = "spearman", exact = FALSE
      )$estimate
      # m3
      predict_day_rank_df[cur + 2, "spearman"] <- cor.test(real_gl_vec,
        rowSums(est_intensity_array_m3[d_test - 15, s, , ]),
        method = "spearman", exact = FALSE
      )$estimate
      # mmhp
      predict_day_rank_df[cur + 3, "spearman"] <- cor.test(real_gl_vec,
        rowSums(est_intensity_array_mmhp[d_test - 15, s, , ]),
        method = "spearman", exact = FALSE
      )$estimate

      # DSNL (not related to observation windows)
      ## updated here also
      dsnl_lambda <- sim_predict_dsnl %>% select(starts_with("lambda_d"))
      dsnl_lambda_arr <- array(as.matrix(dsnl_lambda),
                               dim = c(1000, 6, 12, 12))
      est_intensity_array_dsnl[d_test - 15, s, , ] <- 
        dsnl_lambda_arr[s, d_test - 15, , ]
      predict_day_rank_df[cur + 4, "spearman"] <- cor.test(real_gl_vec,
        rowSums(est_intensity_array_dsnl[d_test - 15, s, , ]),
        method = "spearman", exact = FALSE
      )$estimate

      predict_day_rank_df[c(cur:(cur + 4)), "day"] <- rep(d_test, 5)
      predict_day_rank_df[c(cur:(cur + 4)), "method"] <- c("m1",
                                                           "m2",
                                                           "m3",
                                                           "mmhp",
                                                           "dsnl")
      cur <- cur + 5
    }
    all_cohort_rank_df[cur_all,
                       "spearman"] <- cor.test(real_gl_vec,
                                               rowSums(apply(est_intensity_array_m1[d_test - 15, , , ], c(2, 3), median)), exact = FALSE)$estimate
    all_cohort_rank_df[cur_all + 1,
                       "spearman"] <- cor.test(real_gl_vec,
                                               rowSums(apply(est_intensity_array_m2[d_test - 15, , , ], c(2, 3), median)), exact = FALSE)$estimate
    all_cohort_rank_df[cur_all + 2,
                       "spearman"] <- cor.test(real_gl_vec,
                                               rowSums(apply(est_intensity_array_m3[d_test - 15, , , ], c(2, 3), median)), exact = FALSE)$estimate
    all_cohort_rank_df[cur_all + 3,
                       "spearman"] <- cor.test(real_gl_vec,
                                               rowSums(apply(est_intensity_array_mmhp[d_test - 15, , , ], c(2, 3), median)), exact = FALSE)$estimate
    all_cohort_rank_df[cur_all + 4,
                       "spearman"] <- cor.test(real_gl_vec,
                                               rowSums(apply(est_intensity_array_dsnl[d_test - 15, , , ], c(2, 3), median)), exact = FALSE)$estimate
    all_cohort_rank_df[c(cur_all:(cur_all + no_method - 1)),
                       "day"] <- rep(d_test, no_method)
    all_cohort_rank_df[c(cur_all:(cur_all + no_method - 1)),
                       "cohort"] <- rep(current_cohort, no_method)
    all_cohort_rank_df[c(cur_all:(cur_all + no_method - 1)),
                       "method"] <- c("m1", "m2", "m3", "mmhp", "dsnl")
    cur_all <- cur_all + no_method
  }
  predict_day_rank_df_lst[[current_cohort]] <- predict_day_rank_df
}


#### save the output ####
save(predict_day_rank_df_lst, all_cohort_rank_df,
  file = paste(data_path, "plot_rank_predict.RData", sep = "")
)
