# library(here)
# library(tidyverse)

plot.s <- 26 # 8 is pretty good
plot.i <- 10
plot.j <- 9

# 29 10 -> 9 also decent
# 26 10 -> 9 also decent
#17 10 -> 9 decent
#1, 10 -> 8 not too bad
#8 10 -> 9 also

### only going up to 100 now
cur_i <- plot.i
cur_j <- plot.j

sim_data_path <- here("output", "revisions", "lapl_check", "sim_lapl", "/")

load(paste0(
  sim_data_path, "/", "sim_model3_fit123_",
  plot.s,
  ".RData"
))
test.mmhp <- sim_model3_data$mmhp_matrix[plot.i, plot.j][[1]]
temp.t <- test.mmhp$tau
current.n <- length(temp.t) - 1
time.segment <- seq(0, tail(temp.t, 1), length.out = 5000)

### Preprocess the model inference result
## model 1
#########

## get f from cmdstan 

f_draws <- sim_model3_stan_sim1 %>% 
  select(starts_with("f"))

model1_par_est <- list(
  lambda0 = mean(sim_model3_stan_sim1$lambda0),
  eta_1 = mean(sim_model3_stan_sim1$eta_1),
  eta_2 = mean(sim_model3_stan_sim1$eta_2),
  eta_3 = mean(sim_model3_stan_sim1$eta_3),
  beta = mean(sim_model3_stan_sim1$beta),
  f = apply(f_draws, 2, mean)
)
lambda.m1 <- uniHawkesIntensityNumeric(
  object = list(
    lambda0 =
      model1_par_est$lambda0,
    alpha =
      model1_fn$alpha.fun(
        model1_par_est$f[plot.i],
        model1_par_est$f[plot.j],
        model1_par_est$eta_1,
        model1_par_est$eta_2,
        model1_par_est$eta_3
      ),
    beta = model1_par_est$beta
  ),
  events = temp.t,
  time.vec = time.segment
)


## model2 ####
gamma_draws <- sim_model3_stan_sim2 %>% 
  select(starts_with("gamma"))

zeta_draws <- sim_model3_stan_sim2 %>% 
  select(starts_with("zeta"))

f_draws <- sim_model3_stan_sim2 %>% 
  select(starts_with("f"))

model2_par_est <- list(
  gamma = apply(gamma_draws, 2, mean),
  zeta = apply(zeta_draws, 2, mean),
  eta_1 = mean(sim_model3_stan_sim2$eta_1),
  eta_2 = mean(sim_model3_stan_sim2$eta_2),
  eta_3 = mean(sim_model3_stan_sim2$eta_3),
  beta = mean(sim_model3_stan_sim2$beta),
  f = apply(f_draws, 2, mean)
)

lambda.m2 <- uniHawkesIntensityNumeric(
  object = list(
    lambda0 = model2_par_est$gamma[plot.i] +
      model2_par_est$zeta[plot.j],
    alpha = model1_fn$alpha.fun(
      model2_par_est$f[plot.i],
      model2_par_est$f[plot.j],
      model2_par_est$eta_1,
      model2_par_est$eta_2,
      model2_par_est$eta_3
    ),
    beta = model2_par_est$beta
  ),
  events = temp.t,
  time.vec = time.segment
)

## Model3 ####

lam0_draws <- sim_model3_stan_sim3 %>% 
  select(starts_with("lambda0"))

lam1_draws <- sim_model3_stan_sim3 %>% 
  select(starts_with("lambda1"))

gamma_draws <- sim_model3_stan_sim3 %>% 
  select(starts_with("gamma"))

zeta_draws <- sim_model3_stan_sim3 %>% 
  select(starts_with("zeta"))

f_draws <- sim_model3_stan_sim2 %>% 
  select(starts_with("f"))


lam0_vec <- apply(lam0_draws, 2, mean)
lam1_vec <- apply(lam1_draws, 2, mean)
lam0_par_est <- matrix(0,
                       nrow = length(object_par$f_vec_1),
                       ncol = length(object_par$f_vec_1)
)
lam1_par_est <- matrix(0,
                       nrow = length(object_par$f_vec_1),
                       ncol = length(object_par$f_vec_1)
)

clean_sim_data <- cleanSimulationData(
  raw_data = sim_model3_data,
  cut_off = cut_off,
  N = length(object_par$f_vec_1)
)

for (i in seq_along(lam0_vec)) {
  row_id <- clean_sim_data$I_fit[i]
  col_id <- clean_sim_data$J_fit[i]
  lam0_par_est[row_id, col_id] <- lam0_vec[i]
  lam1_par_est[row_id, col_id] <- lam1_vec[i]
}
###
### model 3 ####
mmhp_par_est <- list(
  lambda0 = lam0_par_est,
  lambda1 = lam1_par_est,
  eta_1 = mean(sim_model3_stan_sim3$eta_1),
  eta_2 = mean(sim_model3_stan_sim3$eta_2),
  eta_3 = mean(sim_model3_stan_sim3$eta_3),
  beta = mean(sim_model3_stan_sim3$beta),
  f = apply(f_draws, 2, mean)
)

object_hat <- list(
  lambda0 = mmhp_par_est$lambda0[cur_i, cur_j],
  lambda1 = mmhp_par_est$lambda1[cur_i, cur_j],
  alpha = model3_fn$alpha.fun(
    mmhp_par_est$f[cur_i],
    mmhp_par_est$f[cur_j],
    mmhp_par_est$eta_1,
    mmhp_par_est$eta_2
  ),
  beta = mmhp_par_est$beta,
  q1 = model3_fn$q1.fun(
    mmhp_par_est$f[cur_i],
    mmhp_par_est$f[cur_j],
    mmhp_par_est$eta_3
  ),
  q2 = model3_fn$q0.fun(
    mmhp_par_est$f[cur_i],
    mmhp_par_est$f[cur_j],
    mmhp_par_est$eta_3
  )
)

## then need to load in the state estimate here
# load(paste0(here("output", "revisions", "sim_m3"),
#   "/", "fit123_state_est_",
#   plot.s,
#   ".RData"
# ))

load(paste0(
  sim_data_path, "/", "fit123_state_est_",
  plot.s,
  ".RData"
))

state.est.latent.mmhp <- interpolate_state_est_lst[plot.i, plot.j][[1]]
state.est.latent.mmhp.new <- fixInterpolateState(state.est.latent.mmhp,
                                                 termination = 100
)

step.fun.est <- stepfun(
  state.est.latent.mmhp$x.hat,
  2 - state.est.latent.mmhp$z.hat
)

lambda.m3 <- mmhpIntensityNumeric(
  params = object_hat,
  t = temp.t[-1],
  time.vec = time.segment,
  latent.vec = step.fun.est(time.segment)
)

object_true <- lapply(object_matrix, function(x) x[plot.i, plot.j])
names(object_true) <- c("lambda0", "lambda1", "alpha", "beta", "q1", "q2")

lambda.true <- mmhpTrueIntensityNumeric(
  params = object_true,
  t = temp.t,
  latent = list(
    x = fixStateTransition(test.mmhp)$x,
    z = fixStateTransition(test.mmhp)$z
  ),
  time.vec = time.segment
)

mmhp_lty <- 1

y.ub <- c(5, 5, 5) # c(130,72,31)
x_events <- 90 # 135#150
my.xlim <- 100 # 160#200#130
my.long.x <- 100 # 205#134.5
legend.x <- c(50, 50, 50) # c(135,135,135)
legend.y <- # c(34,32,32)
legend.y <- c(4, 4, 4) # c(15,12,15)#c(62,37,18)#c(22,22,22)
legend.cex <- 1.5 # 3
line.alpha <- c(1, 0.5)
line.wdth <- c(1.5, 2)
negative.col <- add.alpha("lightskyblue")
positive.col <- add.alpha("tomato")


model_colors <- col_df %>% filter(method %in% c(
  "C-HP",
  "C-DCHP",
  "C-MMHP"
))

model_colors <- model_colors %>% pull(cols)

# model1
############
true.object <- object_true


layout(matrix(c(1:6), 6, 1), heights = rep(c(7.5, 2.5), 3))
par(
  mar = c(0.2, 0.3, 0, 0.1), oma = c(3, 3.5, 4, 2),
  tcl = 0.2, mgp = c(0.5, 0, 0), xpd = TRUE
)

mmhp_est <- fixStateTransition(test.mmhp)

drawUniMMHPIntensityPaper(true.object,
                          simulation = mmhp_est,
                          yupper = y.ub[1],
                          color = add.alpha("black", alpha = line.alpha[1]),
                          line.width = line.wdth[1],
                          y.ratio = -0.05, min.y = -0.5,
                          min.x = 0, max.x = my.xlim,
                          title_name = paste(plot.i, "->", plot.j),
                          title.cex = 2,
                          box.type = "n",
                          line_type = mmhp_lty
)
drawHawkesIntensityPaper(
  lambda0 = model1_par_est$lambda0,
  alpha = model1_fn$alpha.fun(
    model1_par_est$f[plot.i],
    model1_par_est$f[plot.j],
    model1_par_est$eta_1,
    model1_par_est$eta_2,
    model1_par_est$eta_3
  ),
  beta = model1_par_est$beta,
  events = test.mmhp$tau,
  color = add.alpha(model_colors[1],
                    alpha = line.alpha[2]
  ),
  line.width = line.wdth[2]
)

# draw box
axis(1,
     at = c(-2, my.long.x),
     col = "gray35",
     line = 0.5,
     tick = T,
     labels = rep("", 2),
     lwd = 0.5,
     lwd.ticks = 0,
     lty = 3
)


## add legend
# changed 34 to 23, 25 to 18, -20 to -10
# legend(10,legend.y[1],c("State 1/0 events","State change point"),
#        col = c(NA,"red"),
#        y.intersp=0.85,x.intersp=-0.1, bty = "n",
#        pch = c(NA,4), pt.cex = c(NA,2), cex=legend.cex,
#        lty = c(NA,NA), lwd=c(NA,4))
# points(29,84,pch=1,cex=2,col="blue")
# points(26,84,pch=16,cex=2,col="blue")

#### more legend
legend(legend.x[1] - 2, legend.y[1] + 1, c("True"),
       col = c("black"),
       y.intersp = 0.88, x.intersp = 0.65, bty = "n",
       lty = c(1), lwd = c(5), cex = legend.cex, seg.len = 1
)
legend(legend.x[1] + 22, legend.y[1] + 1, c("C-HP"),
       col = c(model_colors[1]),
       y.intersp = 0.88, x.intersp = 0.65, bty = "n",
       lty = c(1), lwd = c(5), cex = legend.cex, seg.len = 1
)

## delta lambda
par(mar = c(0.1, 0.3, 0, 0.1))
plot(0, 0,
     xlim = c(0, my.xlim), 
     ylim = c(-1.5, 1.5),
     type = "n",
     bty = "n",
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n",
     axes = FALSE
)
delta.lambda <- lambda.m1$lambda.t - lambda.true$lambda.t
m1_diff <- delta.lambda
delta.positive <- delta.lambda
delta.positive[delta.positive < 0] <- 0
delta.negative <- delta.lambda
delta.negative[delta.negative > 0] <- 0
### hiding these for now ###
polygon(c(0, lambda.true$time.vec, x_events), # 200
        c(0, delta.negative, 0),
        col = negative.col, border = NA
)
polygon(c(0, lambda.true$time.vec, x_events), # 200
        c(0, delta.positive, 0),
        col = positive.col, border = NA
)

axis(1,
     at = c(-1, my.long.x), 
     col = "gray35",
     line = 0, tick = T, 
     labels = rep("", 2),
     lwd = 0.7,
     lwd.ticks = 0
)

# model2
############
drawUniMMHPIntensityPaper(true.object,
                          simulation = mmhp_est,
                          ### these points aren't consistent?
                          yupper = y.ub[2],
                          color = add.alpha("black",
                                            alpha = line.alpha[1]),
                          line.width = line.wdth[1],
                          y.ratio = -0.05,
                          min.y = -0.5,
                          min.x = 0,
                          max.x = my.xlim,
                          box.type = "n",
                          line_type = mmhp_lty
)
drawHawkesIntensityPaper(
  lambda0 = model2_par_est$gamma[plot.i] +
    model2_par_est$zeta[plot.j],
  alpha = model1_fn$alpha.fun(
    model2_par_est$f[plot.i],
    model2_par_est$f[plot.j],
    model2_par_est$eta_1,
    model2_par_est$eta_2,
    model2_par_est$eta_3
  ),
  beta = model2_par_est$beta,
  events = test.mmhp$tau,
  color = add.alpha(model_colors[2],
                    alpha = line.alpha[2]
  ),
  line.width = line.wdth[2]
)

### used legend
legend(legend.x[2] - 2, legend.y[2] + 1, c("True"),
       col = c("black"),
       y.intersp = 0.88, x.intersp = 0.65, bty = "n",
       lty = c(1), lwd = c(5), cex = legend.cex, seg.len = 1
)
legend(legend.x[2] + 22, legend.y[2] + 1, "C-DCHP",
       col = model_colors[2],
       y.intersp = 0.88,
       x.intersp = 0.65,
       bty = "n",
       lty = 1,
       lwd = 5,
       cex = legend.cex,
       seg.len = 1
)

axis(1,
     at = c(-2, my.long.x),
     col = "gray35",
     line = 0.5,
     tick = T,
     labels = rep("", 2),
     lwd = 0.5,
     lwd.ticks = 0,
     lty = 3
)


## location of legend
legend(legend.x[2] - 35,
       legend.y[2] + 1.75, c("State 0/1 events", "State change point"),
       col = c(NA, "red"),
       y.intersp = 0.35,
       x.intersp = -0,
       bty = "n",
       pch = c(NA, 4), pt.cex = c(NA, 2), cex = legend.cex,
       lty = c(NA, NA), lwd = c(NA, 2)
)
points(legend.x[2] - 31, legend.y[2] + 0.15, 
       pch = 1, cex = 2, col = "blue")
points(legend.x[2] - 30, legend.y[2] + 0.15,
       pch = 16, cex = 2, col = "blue")

## delta lambda
par(mar = c(0.1, 0.3, 0, 0.1))
delta.lambda <- lambda.m2$lambda.t - lambda.true$lambda.t

plot(lambda.true$time.vec, delta.lambda,
     xlim = c(0, my.xlim),
     ylim = c(-1.5, 1.5), type = "n",
     bty = "n",
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n",
     axes = FALSE
)
delta.lambda <- lambda.m2$lambda.t - lambda.true$lambda.t
m2_diff <- delta.lambda
delta.positive <- delta.lambda
delta.positive[delta.positive < 0] <- 0
delta.negative <- delta.lambda
delta.negative[delta.negative > 0] <- 0
polygon(c(0, lambda.true$time.vec, x_events), # 200
        c(0, delta.negative, 0),
        col = negative.col, border = NA
)
polygon(c(0, lambda.true$time.vec, x_events), # 200
        c(0, delta.positive, 0),
        col = positive.col, border = NA
)
# segments(0,0,my.xlim,0,col="gray45")
## draw box
axis(1,
     at = c(-2, my.long.x), col = "gray35",
     line = 0, tick = T, labels = rep("", 2), lwd = 0.7, lwd.ticks = 0
)


# model3
########
par(mar = c(0.1, 0.3, 1, 0.1))
drawUniMMHPIntensityPaper(true.object,
                          simulation = mmhp_est,
                          yupper = y.ub[3],
                          color = add.alpha("black", alpha = line.alpha[1]),
                          line.width = line.wdth[1],
                          y.ratio = -0.05, min.y = -0.5,
                          min.x = 0, max.x = my.xlim,
                          box.type = "n", line_type = mmhp_lty
)
drawUniMMHPIntensityPaper(object_hat,
                          simulation = list(
                            x = state.est.latent.mmhp.new$x.hat,
                            z = state.est.latent.mmhp.new$z.hat,
                            tau = test.mmhp$tau, zt = test.mmhp$zt,
                            lambda.max = test.mmhp$lambda.max
                          ),
                          yupper = y.ub, add = TRUE,
                          color = add.alpha(model_colors[3],
                                            alpha = line.alpha[2]
                          ),
                          line.width = line.wdth[2]
)
legend(legend.x[3] - 2, legend.y[3] + 1.75, c("True"),
       col = c("black"),
       y.intersp = 0.88, x.intersp = 0.65, bty = "n",
       lty = c(1), lwd = c(5), cex = legend.cex, seg.len = 1
)


legend(legend.x[3] + 22, legend.y[3] + 1.75, "C-MMHP",
       col = model_colors[3],
       y.intersp = 0.88,
       x.intersp = .25,
       bty = "n",
       lty = 1,
       lwd = 5,
       cex = legend.cex,
       seg.len = 1
)
## draw box
axis(1,
     at = c(-2, my.long.x),
     col = "gray35",
     line = 0.5,
     tick = T,
     labels = rep("", 2), 
     lwd = 0.5,
     lwd.ticks = 0, lty = 3
)

legend(legend.x[3] - 45, legend.y[3] + 2.5,
       c("Overestimation", "Underestimation"),
       col = c(positive.col, negative.col),
       y.intersp = 0.25, x.intersp = 0.5, bty = "n",
       lty = c(1, 1), lwd = c(10, 10), cex = legend.cex, seg.len = 0.8
)

## delta lambda
par(mar = c(0.1, 0.3, 0, 0.1))
plot(0, 0,
     xlim = c(0, my.xlim), 
     ylim = c(-1.5, 1.5), type = "n",
     bty = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     axes = FALSE
)
delta.lambda <- lambda.m3 - lambda.true$lambda.t
m3_diff <- delta.lambda
delta.positive <- delta.lambda
delta.positive[delta.positive < 0] <- 0
delta.negative <- delta.lambda
delta.negative[delta.negative > 0] <- 0
### hiding these for now
polygon(c(0, lambda.true$time.vec, x_events), # 200
        c(0, delta.negative, 0),
        col = negative.col, border = NA
)
polygon(c(0, lambda.true$time.vec, x_events), # 200
        c(0, delta.positive, 0),
        col = positive.col, border = NA
)
# segments(0,0,my.xlim,0,col="gray45")
## draw box
axis(1,
     at = c(-2, my.long.x), col = "gray35",
     line = 0, tick = T, labels = rep("", 2), lwd = 0.7, lwd.ticks = 0
)


## top
mtext(text = "(b)", side = 3, line = 0.25,
      outer = TRUE, cex = 2.5, font = 2)
## bottom
mtext(text = "Time", side = 1, line = 1.5,
      outer = TRUE, cex = 2.25)
## left
mtext(text = "Intensity", side = 2,
      line = 0.6, outer = TRUE, cex = 2.25)


summary(abs(m1_diff))
summary(abs(m2_diff))
summary(abs(m3_diff))


##### examining the problems with the plot some more ####

#### compare inferred base rates to truth ###
# sim_model3_stan_sim3 %>% summarise_draws()
# 
# sim_model3_stan_sim3 %>% select(starts_with("lambda0")) %>% 
#   apply(., 2, mean)
# 
# true.object
# 
# object_hat
# 
# ##
# mmhp_est
# state.est.latent.mmhp.new
