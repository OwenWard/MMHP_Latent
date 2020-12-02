y.ub <- c(15,15,15)#c(130,72,31)
x_events <- 190#135#150
layout(matrix(c(1:6), 6, 1), heights=rep(c(7.5,1),3))
my.xlim <- 200#160#200#130
my.long.x <- 205#205#134.5
legend.x <- c(150,150,150)#c(135,135,135)
legend.y <- #c(34,32,32)
legend.y <- c(12,12,12)#c(15,12,15)#c(62,37,18)#c(22,22,22)
legend.cex <- 1.5 #3
line.alpha <- c(0.9,0.6)
line.wdth <- c(2,2.5)
negative.col <- add.alpha('lightskyblue')
positive.col <- add.alpha('tomato')
#negative.col <- viridis::viridis(5,alpha = 0.6)[1]
#positive.col <- viridis::viridis(5,alpha = 0.6)[4]
par(mar = c(0,0.3,0,0.1), oma = c(3,3.5,4,0),
    tcl=0.2,mgp=c(0.5,0,0), xpd=TRUE)

# model1
############
true.object <- object_true

drawUniMMHPIntensityPaper(true.object,simulation=fixStateTransition(test.mmhp),
                          yupper=y.ub[1], 
                          color=add.alpha("black", alpha=line.alpha[1]),
                          line.width=line.wdth[1],
                          y.ratio=-0.05, min.y=-2, min.x=0, max.x=my.xlim,
                          title_name=paste(plot.i,"->",plot.j),
                          title.cex = 2,
                          box.type="n")
drawHawkesIntensityPaper(lambda0=model1_par_est$lambda0,
                         alpha=model1_fn$alpha.fun(model1_par_est$f[plot.i],
                                                   model1_par_est$f[plot.j],
                                                   model1_par_est$eta_1,
                                                   model1_par_est$eta_2,
                                                   model1_par_est$eta_3),
                         beta=model1_par_est$beta,
                         events=test.mmhp$tau,
                         color=add.alpha(model_colors[1],
                                         alpha=line.alpha[2]),
                         line.width=line.wdth[2])

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
legend(legend.x[1]+22,legend.y[1],c("C-HP"),
       col = c(model_colors[1]),
       y.intersp=0.88,x.intersp=0.15,bty = "n",
       lty = c(1), lwd = c(5), cex=legend.cex,seg.len=1)

## delta lambda
par(mar = c(0,0.3,0,0.1))
plot(0,0,xlim=c(0,my.xlim), ylim=c(-2,4), type="n",
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

legend(legend.x[2]-65,
       legend.y[2]+2,c("State 0/1 events","State change point"),
       col = c(NA,"red"),
       y.intersp=0.85,x.intersp=-0.1, bty = "n",
       pch = c(NA,4), pt.cex = c(NA,2), cex=legend.cex,
       lty = c(NA,NA), lwd=c(NA,4))
points(90,legend.y[2] - 1.5, pch = 1, cex = 2, col = "blue")
points(93,legend.y[2] - 1.5, pch = 16, cex = 2, col = "blue")

## delta lambda
par(mar = c(0,0.3,0,0.1))
delta.lambda <- lambda.m2$lambda.t-lambda.true$lambda.t

plot(lambda.true$time.vec,delta.lambda,xlim=c(0,my.xlim), 
     ylim=c(-2,4), type="n",
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
       y.intersp=0.6,x.intersp=0.1,bty = "n",
       lty =1, lwd = 5,cex=legend.cex,seg.len=1)
## draw box
axis(1,at=c(-2,my.long.x),col="gray35",
     line=0.5,tick=T,labels=rep("",2),lwd=0.5,lwd.ticks=0,lty=3)

legend(legend.x[3]-65,legend.y[3]+2, c("Overestimation","Underestimation"),
       col = c(positive.col,negative.col),
       y.intersp=0.6,x.intersp=0.5,bty = "n",
       lty = c(1,1), lwd = c(10,10), cex=legend.cex,seg.len=0.8)

## delta lambda
par(mar = c(0,0.3,0,0.1))
plot(0,0,xlim=c(0,my.xlim), ylim=c(-2,4), type="n",
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
mtext(text="Time",side=1,line=2,outer=TRUE,cex=2.25)
## left
mtext(text="Intensity",side=2,line=0.6,outer=TRUE,cex=2.25)
