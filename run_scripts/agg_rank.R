# this file refits the corrected aggregate ranking model
# with unconstrained parameters b and c



### code ###
## run this section if running on the cluster
source("/moto/stats/users/ogw2103/Code/MMHP_Latent/run_scripts/cluster_setup.R")
# ### set cohort_id based on job num
jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
jobid <- as.numeric(jobid)
cohort_id <- jobid
# cohort_id <- 1
#####

save_data_path <- "output/revisions/lapl_check/"


library(cmdstanr)
options(mc.cores = parallel::detectCores())


library(compete)
# library(RColorBrewer)
# library(Hmisc)
# library(wCorr)
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





#### then load in the data here ###

full_data <- readRDS("data/mice.RData")
# A=c9, B=c10, C=c12, D=c15, E=c16, F=c17, G=c18, H=c37, I=c38. J=c45
cohort_names <- paste("cohort", c(9, 10, 12, 15, 16, 17, 18, 37, 38, 45), sep = "")
cohort_short_names <- paste("C", c(9, 10, 12, 15, 16, 17, 18, 37, 38, 45), sep = "")
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
  naive_rank_10[[current_cohort]] <- naiveRankHierarchy(full_data[[cohort_names[current_cohort]]])
  expert_rank_10[[current_cohort]] <- expertRankHierarchy(full_data[[cohort_names[current_cohort]]])
}




current_cohort <- cohort_id
print(paste("Cohort", current_cohort))

#### fit the stan model ####
clean_data <- cleanData(
  raw_data = full_data[[cohort_names[current_cohort]]],
  cut_off = 1, N = 12
)
# changed the cut off here from 0 to 1
agg_model <- cmdstan_model("lib/latent_rank_agg.stan")
fit_agg_rank <- agg_model$sample(
  data = list(n_matrix = clean_data$N_count),
  iter_sampling = 1000,
  chains = 4,
  thin = 4,
  adapt_delta = 0.99
)

# sim_agg_rank <- rstan::extract(fit_agg_rank)
agg_fit <- fit_agg_rank$draws()
sim_agg_rank <- posterior::as_draws_df(agg_fit)

dir.create(paste(save_data_path,
  cohort_names[current_cohort],
  sep = ""
),
recursive = TRUE, showWarnings = FALSE
)
save(fit_agg_rank, sim_agg_rank,
  file = paste(save_data_path, cohort_names[current_cohort],
    "/agg_rank_stan_result_",
    cohort_names[current_cohort],
    ".RData",
    sep = ""
  )
)

#### predictions for this model ####

no_days <- length(unique(full_data[[cohort_names[current_cohort]]]$day))
N_count_train <- cleanDataForTraining(
  all_raw_data =
    full_data[[cohort_names[current_cohort]]],
  cut_off = cut_off,
  N = 12,
  train_period = 15
)$N_count

fit_predict_agg_rank <- agg_model$sample(
  data = list(n_matrix = N_count_train),
  iter_sampling = 1000,
  chains = 4,
  adapt_delta = 0.99
)

# sim_predict_agg_rank <- rstan::extract(fit_predict_agg_rank)
agg_fit <- fit_predict_agg_rank$draws()
sim_predict_agg_rank <- posterior::as_draws_df(agg_fit)

dir.create(paste(save_data_path,
                 cohort_names[current_cohort],
                 sep = ""),
  recursive = TRUE,
  showWarnings = FALSE
)

save(fit_predict_agg_rank, sim_predict_agg_rank,
  file = paste(save_data_path, cohort_names[current_cohort],
    "/predict_agg_rank_stan_result_", cohort_names[current_cohort],
    ".RData",
    sep = ""
  )
)
