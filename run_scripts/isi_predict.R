### Predictions at a cohort level ####
### Here we extract the first 15 days of data for each cohort,
### Determine the rankings and use these to compute the I&SI score,
### based on this ranking, for the remaining events

### code ###
## run this if running on the cluster
# source("/rigel/stats/users/ogw2103/code/MMHP/MMHP_Latent/run_scripts/cluster_setup.R")
# ### set cohort_id based on job num
# jobid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# jobid <- as.numeric(jobid)
# cohort_id <- jobid
#cohort_id <- 1
####
data_path <- "output/"


library(rstan)
library(tidyverse)
options(mc.cores = parallel::detectCores())


library(compete)
library(PlayerRatings)
library(here)
#library(RColorBrewer)
#library(Hmisc)
#library(wCorr)
#library(tidyverse)
library(R.utils)
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






#### then load in the data ####

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


mini_wl_matrix <- function(df, num_mice = 12) {
  # to avoid issues with reordering the data
  df <- as.data.frame(df)
  mylevs <- 1:num_mice
  df[, 1] <- factor(df[, 1], levels = mylevs)
  df[, 2] <- factor(df[, 2], levels = mylevs)
  df$result <- 1
  m1 = stats::xtabs(result ~ ., data = df)
  m1
}


compute_isi_score <- function(win_loss_matrix, curr_ranking) {
  # this takes ranking from highest to lowest
  wl <- win_loss_matrix[curr_ranking,curr_ranking]
  n <- as.numeric(dim(wl)[1])
  m <- matrix(0, n, n)
  siind <- matrix(0, n, n)
  I = 0
  SI = 0
  for (j in 2:n) {
    for (i in 1:(j - 1)) {
      siind[j, i] = j - i
      siind[i, j] = j - i
      if (wl[j, i] == 0 & wl[i, j] == 0) {
        m[j, i] = 0
        m[i, j] = 0
      }
      else {
        if (wl[j, i] > wl[i, j]) {
          m[j, i] = 1
          m[i, j] = -1
        }
        if (wl[j, i] < wl[i, j]) {
          m[j, i] = -1
          m[i, j] = 1
        }
        if (wl[j, i] == wl[i, j]) {
          m[j, i] = 0.5
          m[i, j] = 0.5
        }
      }
    }
  }
  I = sum(lower.tri(m) * m == 1)
  SI = sum((lower.tri(m) * m == 1) * siind)
  
  list(I_score = I, SI_score = SI)
}

# Define the cohorts will be fitted
fit_cohorts <- c(1:10)
naive_rank_10 <- list()
expert_rank_10 <- list()
for(current_cohort in fit_cohorts){
  naive_rank_10[[current_cohort]] <- naiveRankHierarchy(full_data[[cohort_names[current_cohort]]])
  expert_rank_10[[current_cohort]] <- expertRankHierarchy(full_data[[cohort_names[current_cohort]]])
}

for(current_cohort in 1:10) {
  # current_cohort <- cohort_id
  print(paste("Cohort",current_cohort))
  
  training_day <- 15
  ## if I change this then I need to refit
  
  print(current_cohort)
  clean_data <- cleanData(full_data[[cohort_names[current_cohort]]])
  return_df <- cleanObservationPeriod(current_cohort, clean_data)
  
  training_data <- return_df %>%
    filter(day <= training_day)
  
  test_data <- return_df %>%
    filter(day > training_day)
  
  
  #### Fit M1 to half data####
  load(paste(data_path,cohort_names[current_cohort],
             "/cohort_hp_predict_stan_result_",cohort_names[current_cohort],
             ".RData",sep=''))
  
  # highest to lowest
  
  rank_m1 <- rev(order(apply(sim_cohort_hp$f, 2, mean)))
  
  #### Predict M2 ####
  load(paste(data_path,cohort_names[current_cohort],
             "/cohort_dchp_predict_stan_result_",cohort_names[current_cohort],
             ".RData",sep=''))
  # highest to lowest
  rank_m2 <- rev(order(apply(sim_cohort_dchp$f, 2, mean)))
  
  #### Predict M3 ####
  load(paste(data_path,cohort_names[current_cohort],
             "/cohort_mmhp_predict_stan_result_",cohort_names[current_cohort],
             ".RData",sep=''))
  
  # highest to lowest
  rank_m3 <- rev(order(apply(sim_cohort_mmhp$f, 2, mean)))
  
  #### Predict Agg Ranking ####
  load(paste(data_path,cohort_names[current_cohort],
             "/predict_agg_rank_stan_result_",cohort_names[current_cohort],
             ".RData",sep=''))
  
  # highest to lowest
  rank_agg <- rev(order(apply(sim_predict_agg_rank$x, 2, mean)))
  
  #### ISI Predict ####
  train_wl <- mini_wl_matrix(df = cbind(training_data$initiator,
                                        training_data$recipient))
  
  
  isi_train <- isi98(train_wl)
  isi_train_rank <- as.numeric(isi_train$best_order)
  
  #### glicko predict
  glick_df <- training_data %>% 
    mutate(id = row_number(), score = 1) %>%
    select(id, initiator, recipient , score)
  
  gl_train <- my_glicko(glick_df, history=TRUE, cval=2)
  
  gl_train
  
  gl_ranks <- rev(order(gl_train$ratings$Rating))
  
  
  #### save this output ####
  
  ### need to make sure the isi ordering is correct here
  
  est_ranks <- tibble(isi_train = isi_train_rank,
                      agg_train = rank_agg,
                      glicko_train = gl_ranks,
                      m1_train = rank_m1,
                      m2_train = rank_m2,
                      m3_train = rank_m3)
  
  test_wl <- mini_wl_matrix(cbind(test_data$initiator, test_data$recipient))
  best_isi <- isi98(test_wl)
  # this is giving me the ranks correctly but am I mapping them
  # correctly for reordering the wl matrix?
  
  ###
  ## do I want order here?
  ###
  scores <- est_ranks %>% map_dfr(~ compute_isi_score(test_wl, .x))
  
  
  scores$model <- c("isi","agg", "glicko", "m1","m2","m3")
  scores <- scores %>% add_row(tibble_row(I_score = best_isi$I, 
                               SI_score = best_isi$SI,
                               model = "test_isi"))
  
  saveRDS(scores, file = here("output", "rank_sims", "predict_isi",
                              paste0(current_cohort, ".RDS")))
  
}


## then read all these

sim_files <- list.files(here("output","rank_sims","predict_isi"))

all_sims <- list()
for(i in seq_along(sim_files)) {
  single_fit <- readRDS(here("output", "rank_sims", 
                             "predict_isi",
                             sim_files[i]))
  all_sims[[i]] <- single_fit
}

cohort_scores <- all_sims %>%
  map(~as_tibble(.)) %>%
  bind_rows(.id="index")

cohort_scores %>% 
  mutate(ISI = I_score + SI_score) %>%
  ggplot(aes(model,ISI)) + geom_boxplot()

