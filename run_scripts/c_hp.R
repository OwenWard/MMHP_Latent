### fit cohort hawkes process

### code ####
# #run this if running on the cluster
source("/moto/stats/users/ogw2103/Code/MMHP_Latent/run_scripts/cluster_setup.R")
# ### set cohort_id based on job num
jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
cohort_id <- jobid
#cohort_id <- 1
#####

save_data_path <- "output/revisions/aug_run/"

library(cmdstanr)
options(mc.cores = parallel::detectCores())


library(compete)
#library(RColorBrewer)
#library(Hmisc)
#library(wCorr)
#library(tidyverse)
library(dplyr)
#library(R.utils)
#library(fields)



source('lib/naiveRankHierarchy.R')
source('lib/expertRankHierarchy.R')
source('lib/cleanData.R') 
source('lib/prepareDataStan.R') 
source('lib/inferLatentMmhp.R')
source('lib/plotUtil.R') 
source('lib/mmhp.R')
source('lib/uniHawkes.R')
source('lib/simulatePrediction.R')
source('lib/myGlicko.R')
source("https://gist.githubusercontent.com/jalapic/6ca3ece44bdcdc522bb735f183aa0ca0/raw/1a07f469eff08117121b6cbebcd22fb7569e3ee8/compete_extra.R")
source('lib/matrixPlotParameter.R')
source('lib/residualStructureScore.R')





#### then load in the stan fit here ###

full_data <- readRDS("data/mice.RData")
# A=c9, B=c10, C=c12, D=c15, E=c16, F=c17, G=c18, H=c37, I=c38. J=c45
cohort_names <- paste("cohort",c(9,10,12,15,16,17,18,37,38,45),sep='')
cohort_short_names <- paste("C",c(9,10,12,15,16,17,18,37,38,45),sep='')
cut_off <- 3
mice_number <- 12


model1_fn <- list(alpha.fun = function(x,y,eta1,eta2,eta3){return(eta1*x*y*exp(-eta2*abs(x-y))/(1+exp(-eta3*(x-y))))})

model3_fn <- list(alpha.fun = function(x,y,eta1,eta2){return(eta1*x*y*exp(-eta2*abs(x-y)))},
                  q1.fun = function(x,y,eta3){return(exp(-eta3*x))},
                  q0.fun = function(x,y,eta3){return(exp(-eta3*y))})

# Define the cohorts will be fitted
fit_cohorts <- c(1:10)
naive_rank_10 <- list()
expert_rank_10 <- list()
for(current_cohort in fit_cohorts){
  naive_rank_10[[current_cohort]] <- naiveRankHierarchy(full_data[[cohort_names[current_cohort]]])
  expert_rank_10[[current_cohort]] <- expertRankHierarchy(full_data[[cohort_names[current_cohort]]])
}




current_cohort <- cohort_id
print(paste("Cohort",current_cohort))


# then fit C-HP here for each cohort
model1 <- cmdstan_model("lib/model1.stan")

#### fit stan model ####
print(paste("Cohort",current_cohort))
stan_input_lst <- prepareDataStan(current_cohort)
stan_input_lst$alpha_id <- expert_rank_10[[current_cohort]][1]
stan_input_lst$delta_1 <- rep(0.5,stan_input_lst$N_til)




fit_cohort_hp <- model1$sample(data = stan_input_lst,
                               iter_warmup = 1000,
                               iter_sampling = 1000,
                               chains = 4,
                               thin = 4,
                               adapt_delta=0.98,
                               refresh = 200)

chp_fit <- fit_cohort_hp$draws()
sim_cohort_hp <- posterior::as_draws_df(chp_fit)

dir.create(paste(save_data_path,
                 cohort_names[current_cohort],sep=''),
           recursive = TRUE,
           showWarnings = FALSE)
save(sim_cohort_hp, fit_cohort_hp,
     file = paste(save_data_path,
                  cohort_names[current_cohort],
                  "/cohort_hp_stan_result_",
                  cohort_names[current_cohort],
                  ".RData", sep = ''))


### Predictions for this stan model ####
print(current_cohort)
stan_train_input_lst <- prepareDataStanTrain(current_cohort)
stan_train_input_lst$alpha_id <- expert_rank_10[[current_cohort]][1]
stan_train_input_lst$delta_1 <- rep(0.5, stan_train_input_lst$N_til)

fit_cohort_hp <- model1$sample(data = stan_train_input_lst,
                               iter_warmup = 1000, 
                               iter_sampling = 1000,
                               chains = 4,
                               adapt_delta=0.99)



chp_fit <- fit_cohort_hp$draws()
sim_cohort_hp <- posterior::as_draws_df(chp_fit)

# sim_cohort_hp <- rstan::extract(fit_cohort_hp)
dir.create(paste(save_data_path,
                 cohort_names[current_cohort],
                 sep = ''),
           recursive = TRUE,
           showWarnings = FALSE)
save(sim_cohort_hp, fit_cohort_hp,
     file = paste(save_data_path,
                  cohort_names[current_cohort],
                  "/cohort_hp_predict_stan_result_",
                  cohort_names[current_cohort],
                  ".RData",
                  sep = ''))

##### pearson residuals for stan fit ####
mice_number <- 12

print(current_cohort)
clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
return_df <- cleanObservationPeriod(current_cohort,
                                    full_data[[cohort_names[current_cohort]]],
                                    clean_data)
unique_pairs_df <- return_df %>% group_by(initiator, recipient) %>%
  summarize(count=n(),
            observe = list(observe.id),
            observe.length = list(observe.time),
            no.events = list(no.events))
unique_observe_win <- unique(return_df[, c("observe.id","observe.time")])


num_winds <- nrow(unique_observe_win)

load(paste(save_data_path,cohort_names[current_cohort],
           "/cohort_hp_stan_result_",
           cohort_names[current_cohort],
           ".RData",
           sep = ''))

f_draws <- sim_cohort_hp %>% select(starts_with("f"))

model1_par_est <- list(
  lambda0 = mean(sim_cohort_hp$lambda0),
  eta_1 = mean(sim_cohort_hp$eta_1),
  eta_2 = mean(sim_cohort_hp$eta_2),
  eta_3 = mean(sim_cohort_hp$eta_3),
  beta = mean(sim_cohort_hp$beta),
  f = apply(f_draws, 2, mean)
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
m1_residual_matrix <- matrix(0, ncol = mice_number, nrow = mice_number)

m1_residual_array <- array(0, dim = c(mice_number, mice_number, num_winds))

for (i in 1:mice_number) {
  for (j in 1:mice_number) {
    pair <- which(unique_pairs_df$initiator == i & unique_pairs_df$recipient == j)
    if (length(pair) > 0 & (i != j)) {
      par_est <- lapply(model1_par_matrix, function(x) x[i, j])
      names(par_est) <- c("lambda0", "alpha", "beta")
      current_window_vec <- unique_pairs_df$observe[[pair]]
      all_residual <- 0
      for (cur in c(1:num_winds)) { ## check length > 2
        # consider all windows which may or may not have events
        if (cur %in% current_window_vec) {
          # cur_win <- current_window_vec[cur]
          cur_win <- cur
          current_event_time <- return_df[return_df$initiator == i &
            return_df$recipient == j &
            return_df$observe.id == cur_win, "event.times"][[1]]
          current_obs_time <- return_df[return_df$initiator == i &
            return_df$recipient == j &
            return_df$observe.id == cur_win, "observe.time"]
          all_residual <- all_residual + uniHawkesPearsonResidual(
            object = par_est,
            events = current_event_time,
            termination = current_obs_time
          )
          m1_residual_array[i, j, cur] <- uniHawkesPearsonResidual(
            object = par_est,
            events = current_event_time,
            termination = current_obs_time
          )
        }
        else {
          # compute the intensity over empty window
          current_obs_time <- unique_observe_win$observe.time[cur]
          # this is returning null
          # these are all the same now
          all_residual <- all_residual + uniHawkesPearsonResidual(
            object = par_est,
            events = NULL,
            termination = current_obs_time
          )
          m1_residual_array[i, j, cur] <- uniHawkesPearsonResidual(
            object = par_est,
            events = NULL,
            termination = current_obs_time
          )
        }
      }
      m1_residual_matrix[i, j] <- all_residual
    }
  }
}

##### save this mmhp_residual_matrix ####
saveRDS(m1_residual_matrix,
  file = paste(save_data_path, cohort_names[current_cohort],
    "/chp_pr_matrix_", cohort_names[current_cohort],
    ".RDS",
    sep = ""
  )
)

saveRDS(m1_residual_array,
  file = paste(save_data_path, cohort_names[current_cohort],
    "/chp_pr_array_", cohort_names[current_cohort],
    ".RDS",
    sep = ""
  )
)
