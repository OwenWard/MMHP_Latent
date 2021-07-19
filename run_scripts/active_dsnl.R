
### code ####
# #run this if running on the cluster
source("/moto/stats/users/ogw2103/Code/MMHP_Latent/run_scripts/cluster_setup.R")
# ### set cohort_id based on job num
jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
cohort_id <- jobid
#cohort_id <- 1
#####

save_data_path <- "output/revisions/lapl_check/"

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


load(paste(data_path, cohort_names[current_cohort],
           "/cmmhp_est_zt_", cohort_names[current_cohort],
           ".RData",
           sep = ""
))
clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
return_df <- cleanObservationPeriod(current_cohort,
                                    full_data[[cohort_names[current_cohort]]],
                                    clean_data)

unique_pairs_df <- return_df %>%
  group_by(initiator, recipient) %>%
  dplyr::summarize(
    count = n(),
    observe = list(observe.id),
    observe.length = list(observe.time),
    no.events = list(no.events)
  )

unique_observe_win <- unique(return_df[,
                                       c("day",
                                         "observe.id",
                                         "observe.time")])

## 5.  find how many 1s for each pair
print(".......state separation plot........")
total_event_array <- array(0, dim = c(
  mice_number, mice_number,
  max(return_df$observe.id)
))
active_event_array <- array(0, dim = c(
  mice_number, mice_number,
  max(return_df$observe.id)
))

for (i in 1:mice_number) {
  for (j in c(1:mice_number)[-i]) {
    pair <- which(unique_pairs_df$initiator == i & 
                    unique_pairs_df$recipient == j)
    if (length(pair) > 0) {
      current_window_vec <- unique_pairs_df$observe[[pair]]
      for (cur_win in current_window_vec) {
        row_indicator <- return_df$initiator == i & 
          return_df$recipient == j & return_df$observe.id == cur_win
        total_event_array[i, j, cur_win] <- length(return_df[
          row_indicator, "event.times"][[1]])
        active_event_array[i, j, cur_win] <- 
          sum(apply(2 - state_array_list[[pair]][[cur_win]], 1, mean) > 0.5)
      }
    }
  }
}

utility_state_day <- array(0, dim = c(21, mice_number, mice_number))
social_state_day <- array(0, dim = c(21, mice_number, mice_number))

for (cur_day in c(1:21)) {
  cur_wins <- which(unique_observe_win$day == cur_day)
  utility_state_day[cur_day, , ] <- apply(
    active_event_array[, , cur_wins],
    c(1, 2), sum
  )
  social_state_day[cur_day, , ] <- apply(
    total_event_array[, , cur_wins],
    c(1, 2), sum
  ) -
    utility_state_day[cur_day, , ]
}

dsnl_poisson <- cmdstan_model("../lib/dsnl_poisson.stan")

## utility dsnl
fit_dsnl_active <- dsnl_poisson$sample(
  data = list(
    day = 21,
    Gt = utility_state_day,
    c = 1,
    sigma = 1
  ),
  iter_sampling = 1000,
  iter_warmup = 1000,
  chains = 4,
  adapt_delta = 0.95
)

rs <- fit_dsnl_active$summary() %>% pull(rhat)
summary(rs)

active_fit <- fit_dsnl_active$draws()
sim_dsnl_active <- posterior::as_draws_df(active_fit)
### output some summaries here 

# print(fit_dsnl_active, pars = c("x"))
# sim_dsnl_active <- rstan::extract(fit_dsnl_active)

## social dsnl
fit_dsnl_inactive <- dsnl_poisson$sample(
                          data = list(
                            day = 21,
                            Gt = social_state_day,
                            c = 1,
                            sigma = 1
                          ),
                          iter_sampling = 1000,
                          chains = 4,
                          adapt_delta = 0.99
)
# print(fit_dsnl_inactive, pars = c("x"))
# sim_dsnl_inactive <- rstan::extract(fit_dsnl_inactive)

rs <- fit_dsnl_inactive$summary() %>% pull(rhat)
summary(rs)


inactive_fit <- fit_dsnl_inactive$draws()
sim_dsnl_inactive <- posterior::as_draws_df(inactive_fit)

save(sim_dsnl_active, fit_dsnl_active, utility_state_day,
     sim_dsnl_inactive, fit_dsnl_inactive, social_state_day,
     file = paste(data_path, cohort_names[current_cohort],
                  "/dsnl_state_separation_stan_result_",
                  cohort_names[current_cohort],
                  ".RData",
                  sep = ""
     )
)
