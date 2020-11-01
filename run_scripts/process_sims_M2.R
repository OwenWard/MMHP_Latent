library(here)
func_dir <- "/lib/"


source(here("lib",'plotUtil.R')) 
source(here("lib",'drawIntensity.R')) 
source(here("lib",'uniHawkes.R')) 
source(here("lib",'mmhp.R')) 


sim_data_path <- "D:/MMHP_Latent_Sims/M2/"  # local location

sim_num  <- 1
n_sim <- 39 # for now

stan_sims_M1 <- list()
stan_sims_M2 <- list()
stan_sims_M3 <- list()

for(sim_num in seq_len(n_sim)) {
  load(paste(sim_data_path, "sim_model2_fit123_", sim_num,
             ".RData", sep=''))
  stan_sims_M1[[sim_num]] <- sim_model2_stan_sim1
  stan_sims_M2[[sim_num]] <- sim_model2_stan_sim2
  stan_sims_M3[[sim_num]] <- sim_model2_stan_sim3
}

# then have all the posterior draws in these 3 list objects

load(paste(sim_data_path, "sim_model2_", 1,
           ".RData", sep = ""))
## to load the true params
# load(paste(sim_data_path, "sim_model1_fit123_", sim_num,
#            ".RData", sep=''))
# load(paste(sim_data_path, "fit123_state_est_M1_", sim_num,
#            ".RData", sep=''))


# n_sim <- 50


#### Plot Latent Ranks ####

#png(paste(plot_path,"sim_latent_recover.png",sep=""), height=500, width=900)
par(mfrow=c(3, length(object_par$f_vec_1)),
    mar=c(0,0,1,0),
    oma=c(6.5,11,3.5,1))

y.ub <- c(12,12,12,12,12)
text.y <- c(rep(8,4),8)
text.x <- c(0.1,0.2,0.4,0.7,0.9)
text.x.adj <- c(0.3,0.4,0.6,0.5,0.8)
model_colors <- c("steelblue", "goldenrod2", "firebrick2")
model_colors <- viridis::viridis(3)

for(m in c(1:3)){
  for(j in c(1:length(object_par$f_vec_1))){
    plot(1,4,xlim=c(0,1),ylim=c(0,y.ub[j]),
         type="n",xlab="",ylab="",xaxt="n",yaxt="n",axes=FALSE)
    
    ## Posterior density
    # current_result <- eval(parse(text=paste("sim_model3_stan_sim",m,sep="")))
    current_result <- eval(parse(text=paste("stan_sims_M",m,sep="")))
    for(i in c(1:n_sim)){
      lines(density(current_result[[i]]$f[1001:2000,j]),
            col=add.alpha(model_colors[m], alpha=0.8),cex=0.6,lwd=0.8)
    }
    
    ## plot true value
    abline(v=object_par$f_vec_1[j],lwd=1.8,lty=2,col="gray30")
    
    ## plot box
    box.wdt <- 1.5
    box.col <- "black"
    box(lty = 'solid',col=box.col,lwd=0.7)
    
    if(m==3){
      text(text.x.adj[j],text.y[j],text.x[j],cex=2.4)
    }
  }
  ## bottom
  mtext(text=expression(f[1],f[2],f[3],f[4],f[5]), side=1, line=2.4, outer=TRUE, 
        at=c(1:length(object_par$f_vec_1))/length(object_par$f_vec_1)-0.1,
        cex=2.4)
  mtext(text="Latent rank variables", side=1, line=5.3, outer=TRUE, cex=2.7,
        at=0.5)
  ## top
  # mtext(text="(a)", side=3, line=0, outer=TRUE, cex=3,
  #       at=0.5, font=2)
  ## left
  mtext(text=rev(c("C-HP  ","C-DCHP","C-MMHP")), side=2, line=0.5, outer=TRUE, 
        at=c(0:2)/3+0.15,
        cex=2, las=2)
}



#### Then recover latent states ####
### this just uses a single sim from above
model1_fn <- list(alpha.fun = function(x,y,eta1,eta2,eta3){return(eta1*x*y*exp(-eta2*abs(x-y))/(1+exp(-eta3*(x-y))))})

model3_fn <- list(alpha.fun = function(x,y,eta1,eta2){return(eta1*x*y*exp(-eta2*abs(x-y)))},
                  q1.fun = function(x,y,eta3){return(exp(-eta3*x))},
                  q0.fun = function(x,y,eta3){return(exp(-eta3*y))})


# then have all the posterior draws in these 3 list objects
sim_num <- 1
load(paste(sim_data_path, "sim_model1_", sim_num,
           ".RData", sep = ""))
## to load the true params
load(paste(sim_data_path, "sim_model1_fit123_", sim_num,
           ".RData", sep=''))
load(paste(sim_data_path, "fit123_state_est_M1_", sim_num,
           ".RData", sep=''))


# plot.s <- 9#7#9
plot.i <- 4
plot.j <- 5

cur_i <- 4
cur_j <- 5

event_ids <- which(sim_model1_data[[1]]$start == cur_i & 
                     sim_model1_data[[1]]$end == cur_j)
temp.t <- sim_model1_data[[1]]$day_hour[event_ids]
current.n <- length(temp.t)
time.segment <- seq(0,tail(temp.t,1),length.out=5000)

### Preprocess the model inference result
## model 1
#########
model1_par_est <- list(lambda0 =
                         mean(sim_model1_stan_sim1$lambda0[1001:2000]),
                       eta_1=mean(sim_model1_stan_sim1$eta_1[1001:2000]),
                       eta_2=mean(sim_model1_stan_sim1$eta_2[1001:2000]),
                       eta_3=mean(sim_model1_stan_sim1$eta_3[1001:2000]),
                       beta=mean(sim_model1_stan_sim1$beta[1001:2000]),
                       f=apply(sim_model1_stan_sim1$f[1001:2000,],2,mean))
lambda.m1<-uniHawkesIntensityNumeric(object=list(lambda0 = 
                                                   model1_par_est$lambda0,
                                                 alpha =
                                                   model1_fn$alpha.fun(model1_par_est$f[plot.i],
                                                                       model1_par_est$f[plot.j],
                                                                       model1_par_est$eta_1,
                                                                       model1_par_est$eta_2,
                                                                       model1_par_est$eta_3),
                                                 beta=model1_par_est$beta),
                                     events=temp.t,
                                     time.vec=time.segment)


## model2
#############
model2_par_est <- list(gamma=apply(sim_model1_stan_sim2$gamma[1001:2000,],
                                   2,mean),
                       zeta = 
                         apply(sim_model1_stan_sim2$zeta[1001:2000,],2,mean),
                       eta_1=mean(sim_model1_stan_sim2$eta_1[1001:2000]),
                       eta_2=mean(sim_model1_stan_sim2$eta_2[1001:2000]),
                       eta_3=mean(sim_model1_stan_sim2$eta_3[1001:2000]),
                       beta=mean(sim_model1_stan_sim2$beta[1001:2000]),
                       f=apply(sim_model1_stan_sim2$f[1001:2000,],2,mean))

lambda.m2 <- uniHawkesIntensityNumeric(object=list(
  lambda0=model2_par_est$gamma[plot.i]+
    model2_par_est$zeta[plot.j],                                        
  alpha=model1_fn$alpha.fun(model2_par_est$f[plot.i],
                            model2_par_est$f[plot.j],
                            model2_par_est$eta_1,
                            model2_par_est$eta_2,
                            model2_par_est$eta_3),
  beta=model2_par_est$beta),
  events=temp.t,
  time.vec=time.segment)

## Model3
#############
mmhp_par_est <- list(lambda0=mean(sim_model1_stan_sim3$lambda0[1001:2000]),
                     lambda1=mean(sim_model1_stan_sim3$lambda1[1001:2000]),
                     eta_1=mean(sim_model1_stan_sim3$eta_1[1001:2000]),
                     eta_2=mean(sim_model1_stan_sim3$eta_2[1001:2000]),
                     eta_3=mean(sim_model1_stan_sim3$eta_3[1001:2000]),
                     beta=mean(sim_model1_stan_sim3$beta[1001:2000]),
                     f=apply(sim_model1_stan_sim3$f[1001:2000,],2,mean))

object_hat <- list(lambda0=mmhp_par_est$lambda0,
                   lambda1=mmhp_par_est$lambda1,
                   alpha=model3_fn$alpha.fun(mmhp_par_est$f[cur_i],
                                             mmhp_par_est$f[cur_j],
                                             mmhp_par_est$eta_1,
                                             mmhp_par_est$eta_2),
                   beta=mmhp_par_est$beta,
                   q1=model3_fn$q1.fun(mmhp_par_est$f[cur_i],
                                       mmhp_par_est$f[cur_j],
                                       mmhp_par_est$eta_3),
                   q2=model3_fn$q0.fun(mmhp_par_est$f[cur_i],
                                       mmhp_par_est$f[cur_j],
                                       mmhp_par_est$eta_3))


state.est.latent.mmhp <- interpolate_state_est_lst[plot.i,plot.j][[1]]
state.est.latent.mmhp.new <- fixInterpolateState(state.est.latent.mmhp,
                                                 termination=200)

step.fun.est <- stepfun(state.est.latent.mmhp$x.hat,
                        2-state.est.latent.mmhp$z.hat)

lambda.m3 <- mmhpIntensityNumeric(params=object_hat,
                                  t=temp.t,
                                  time.vec=time.segment,
                                  latent.vec=step.fun.est(time.segment))

object_true <- lapply(object_matrix,function(x) x[plot.i,plot.j])
names(object_true) <- c("lambda0","lambda1","alpha","beta")

lambda.true <- hawkesIntensityNumeric(params=object_true,
                                      t=temp.t,
                                      # latent=list(x=fixStateTransition(test.mmhp)$x,
                                      #             z=fixStateTransition(test.mmhp)$z),
                                      time.vec=time.segment)

y.ub <- c(20,20,20)#c(130,72,31)
x_events <- 135#150
#model_colors <- c("steelblue", "goldenrod2", "firebrick2")
model_colors <- viridis::viridis(3)
layout(matrix(c(1:6), 6, 1), heights=rep(c(2.5,1),3))
layout(matrix(c(1:6), 6, 1), heights=rep(c(7.5,1),3))
my.xlim <- 140#160#200#130
my.long.x <- 140#205#134.5
legend.x <- c(50,50,50)#c(135,135,135)
legend.y <- #c(34,32,32)
  legend.y <- c(20,12,15)#c(62,37,18)#c(22,22,22)
legend.cex <- 1.5 #3
line.alpha <- c(0.9,0.6)
line.wdth <- c(2,2.5)
negative.col <- add.alpha('lightskyblue',0.6)
positive.col <- add.alpha('tomato',0.6)
#negative.col <- viridis::viridis(5,alpha = 0.6)[1]
#positive.col <- viridis::viridis(5,alpha = 0.6)[4]
par(mar = c(0,0.3,0,0.1), oma = c(3,3.5,4,0),
    tcl=0.2,mgp=c(0.5,0,0), xpd=TRUE)

# model1
############
true.object <- object_true

# drawUniMMHPIntensityPaper(true.object,simulation=fixStateTransition(test.mmhp),
#                           yupper=y.ub[1], 
#                           color=add.alpha("black", alpha=line.alpha[1]),
#                           line.width=line.wdth[1],
#                           y.ratio=-0.05, min.y=-2, min.x=0, max.x=my.xlim,
#                           title_name=paste(plot.i,"->",plot.j),
#                           title.cex = 2,
#                           box.type="n")
# plot.new()
# drawHawkesIntensityPaper(lambda0 = true.object$lambda0,
#                          alpha = true.object$alpha,
#                          beta = true.object$beta, events = temp.t, 
#                          color = add.alpha("black", alpha=line.alpha[1]),
#                          line.width=line.wdth[1])
ppdiag::drawHPIntensity(hp = true.object, events = temp.t)
est_hp <- list(lambda0=model1_par_est$lambda0,
               alpha=model1_fn$alpha.fun(model1_par_est$f[plot.i],
                                         model1_par_est$f[plot.j],
                                         model1_par_est$eta_1,
                                         model1_par_est$eta_2,
                                         model1_par_est$eta_3),
               beta=model1_par_est$beta)
ppdiag::drawHPIntensity(hp = est_hp, events = temp.t,
                        color = add.alpha(model_colors[1],
                                          alpha=line.alpha[2]),
                        add = TRUE)
# drawHawkesIntensityPaper(lambda0=model1_par_est$lambda0,
#                          alpha=model1_fn$alpha.fun(model1_par_est$f[plot.i],
#                                                    model1_par_est$f[plot.j],
#                                                    model1_par_est$eta_1,
#                                                    model1_par_est$eta_2,
#                                                    model1_par_est$eta_3),
#                          beta=model1_par_est$beta,
#                          events=temp.t,
#                          color=add.alpha(model_colors[1],
#                                          alpha=line.alpha[2]),
#                          line.width=line.wdth[2])

## draw box
axis(1,at=c(-2,my.long.x),col="gray35",
     line=0.5,tick=T,labels=rep("",2),lwd=0.5,lwd.ticks=0,lty=3)


## add legend
# changed 34 to 23, 25 to 18, -20 to -10
# legend(10,legend.y[1],c("State 1/0 events","State change point"),
#        col = c(NA,"red"),
#        y.intersp=0.85,x.intersp=-0.1, bty = "n",
#        pch = c(NA,4), pt.cex = c(NA,2), cex=legend.cex,
#        lty = c(NA,NA), lwd=c(NA,4))
# points(29,84,pch=1,cex=2,col="blue")
# points(26,84,pch=16,cex=2,col="blue")
legend(legend.x[1],legend.y[1],c("True"),
       col = c("black"),
       y.intersp=0.88,x.intersp=0.15,bty = "n",
       lty = c(1), lwd = c(5), cex=legend.cex,seg.len=1)
legend(legend.x[1]+18,legend.y[1],c("C-HP"),
       col = c(model_colors[1]),
       y.intersp=0.88,x.intersp=0.15,bty = "n",
       lty = c(1), lwd = c(5), cex=legend.cex,seg.len=1)

## delta lambda
par(mar = c(0,0.3,0,0.1))
plot(0,0,xlim=c(0,my.xlim), ylim=c(-4,10), type="n",
     bty="n", xlab="",ylab="",xaxt="n",yaxt="n",axes=FALSE)
delta.lambda <- lambda.m1$lambda.t-lambda.true$lambda.t
delta.positive <- delta.lambda
delta.positive[delta.positive<0] <- 0
delta.negative <- delta.lambda
delta.negative[delta.negative>0] <- 0
### hiding these for now ###
polygon(c(0,lambda.true$time.vec,x_events),#200
        c(0,delta.negative,0),col=negative.col,border=NA)
polygon(c(0,lambda.true$time.vec,x_events),#200
        c(0,delta.positive,0),col=positive.col,border=NA)

axis(1,at=c(-2,my.long.x),col="gray35",
     line=0,tick=T,labels=rep("",2),lwd=0.7,lwd.ticks=0)


# model2
############
drawUniMMHPIntensityPaper(true.object,simulation=fixStateTransition(test.mmhp),
                          yupper=y.ub[2], 
                          color=add.alpha("black", alpha=line.alpha[1]),
                          line.width=line.wdth[1],
                          y.ratio=-0.05, min.y=-2, min.x=0, max.x=my.xlim, 
                          box.type="n")
drawHawkesIntensityPaper(lambda0=model2_par_est$gamma[plot.i]+
                           model2_par_est$zeta[plot.j],
                         alpha=model1_fn$alpha.fun(model2_par_est$f[plot.i],
                                                   model2_par_est$f[plot.j],
                                                   model2_par_est$eta_1,
                                                   model2_par_est$eta_2,
                                                   model2_par_est$eta_3),
                         beta=model2_par_est$beta,
                         events=test.mmhp$tau,
                         color=add.alpha(model_colors[2],
                                         alpha=line.alpha[2]),
                         line.width=line.wdth[2])
legend(legend.x[2],legend.y[2],"C-DCHP",
       col = model_colors[2],
       y.intersp=0.88,x.intersp=0.1,bty = "n",
       lty = 1, lwd = 5, cex= legend.cex,seg.len=1)

axis(1,at=c(-2,my.long.x),col="gray35",
     line=0.5,tick=T,labels=rep("",2),lwd=0.5,lwd.ticks=0,lty=3)

legend(10,legend.y[2],c("State 0/1 events","State change point"),
       col = c(NA,"red"),
       y.intersp=0.85,x.intersp=-0.1, bty = "n",
       pch = c(NA,4), pt.cex = c(NA,2), cex=legend.cex,
       lty = c(NA,NA), lwd=c(NA,4))
points(14,8,pch=1,cex=2,col="blue")
points(16,8,pch=16,cex=2,col="blue")

## delta lambda
par(mar = c(0,0.3,0,0.1))
delta.lambda <- lambda.m2$lambda.t-lambda.true$lambda.t

plot(lambda.true$time.vec,delta.lambda,xlim=c(0,my.xlim), ylim=c(-3,10), type="n",
     bty="n", xlab="",ylab="",xaxt="n",yaxt="n",axes=FALSE)
delta.lambda <- lambda.m2$lambda.t-lambda.true$lambda.t
delta.positive <- delta.lambda
delta.positive[delta.positive<0] <- 0
delta.negative <- delta.lambda
delta.negative[delta.negative>0] <- 0
polygon(c(0,lambda.true$time.vec,x_events), #200
        c(0,delta.negative,0),col=negative.col,border=NA)
polygon(c(0,lambda.true$time.vec,x_events), #200
        c(0,delta.positive,0),col=positive.col,border=NA)
# segments(0,0,my.xlim,0,col="gray45")
## draw box
axis(1,at=c(-2,my.long.x),col="gray35",
     line=0,tick=T,labels=rep("",2),lwd=0.7,lwd.ticks=0)


# model3
########
par(mar = c(0,0.3,1,0.1))
drawUniMMHPIntensityPaper(true.object,
                          simulation=fixStateTransition(test.mmhp),
                          yupper=y.ub[3],
                          color=add.alpha("black", alpha=line.alpha[1]),
                          line.width=line.wdth[1],
                          y.ratio=-0.05, min.y=-2, min.x=0, max.x=my.xlim, 
                          box.type="n")
drawUniMMHPIntensityPaper(object_hat,
                          simulation = list(
                            x=state.est.latent.mmhp.new$x.hat,                                          
                            z=state.est.latent.mmhp.new$z.hat,
                            tau=test.mmhp$tau,zt=test.mmhp$zt,
                            lambda.max=test.mmhp$lambda.max),
                          yupper=y.ub,add=TRUE,
                          color=add.alpha(model_colors[3],
                                          alpha=line.alpha[2]),
                          line.width=line.wdth[2])
legend(legend.x[3],legend.y[3],"C-MMHP",
       col = model_colors[3],
       y.intersp=0.88,x.intersp=0.1,bty = "n",
       lty =1, lwd = 5,cex=legend.cex,seg.len=1)
## draw box
axis(1,at=c(-2,my.long.x),col="gray35",
     line=0.5,tick=T,labels=rep("",2),lwd=0.5,lwd.ticks=0,lty=3)

legend(10,legend.y[3],c("Overestimation","Underestimation"),
       col = c(positive.col,negative.col),
       y.intersp=0.88,x.intersp=0.5,bty = "n",
       lty = c(1,1), lwd = c(10,10), cex=legend.cex,seg.len=0.8)

## delta lambda
par(mar = c(0,0.3,0,0.1))
plot(0,0,xlim=c(0,my.xlim), ylim=c(-12,8), type="n",
     bty="n", xlab="",ylab="",xaxt="n",yaxt="n",axes=FALSE)
delta.lambda <- lambda.m3-lambda.true$lambda.t
delta.positive <- delta.lambda
delta.positive[delta.positive<0] <- 0
delta.negative <- delta.lambda
delta.negative[delta.negative>0] <- 0
### hiding these for now
polygon(c(0,lambda.true$time.vec,x_events), # 200 
        c(0,delta.negative,0),col=negative.col,border=NA)
polygon(c(0,lambda.true$time.vec,x_events), # 200 
        c(0,delta.positive,0),col=positive.col,border=NA)
# segments(0,0,my.xlim,0,col="gray45")
## draw box
axis(1,at=c(-2,my.long.x),col="gray35",
     line=0,tick=T,labels=rep("",2),lwd=0.7,lwd.ticks=0)


## top
mtext(text="(b)",side=3,line=0.8,outer=TRUE,cex=2.5,font=2)
## bottom
mtext(text="Time",side=1,line=2,outer=TRUE,cex=2.5)
## left
mtext(text="Intensity",side=2,line=0.6,outer=TRUE,cex=2.5)

