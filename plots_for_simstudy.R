# Code to create figures for results from simulation study
# Olivia Bernstein Morgan

library(RColorBrewer)

# round variable and set the numebr of trailing 0's
rnd = function(x, n.dig){format(round(x, digits = n.dig), nsmall = n.dig)}

smooth.plot = function(x,y, n.point = 100, alpha = 0.8){
  # takes x and y and returns points for a smoothed line
  # with n.point number of points
  
  x.smooth = seq(min(x),max(x),length.out = n.point)
  y.smooth = predict(loess(y~x, span = alpha), newdata = x.smooth)
  out = data.frame(x = x.smooth, y = y.smooth)
  return(out)
}


make.plot.smooth.coefX = function(coef.save, 
                                  ylim.percent = NULL,
                                  ylim.absolute = NULL,
                                  title = NULL, 
                                  plot.points = TRUE,
                                  xlab = "Coefficient Magnitude",
                                  ylab = "Absolute Bias (vs SRS)",
                                  alpha1 = .9,
                                  alpha2 = .9,
                                  beta.list, xaxt = NULL){
  # get percent and absolute bias vs SRS
  SRS = matrix(rep(coef.save[,1],6), nrow = length(beta.list),ncol = 6)
  BS = coef.save[,2:7]
  
  # Percent.bias = abs((BS-SRS)/SRS*100)
  Absolute.bias = abs(BS-SRS)
  
  if(is.null(ylim.absolute)){
    ylim.absolute = c(0,max(Absolute.bias))
  }
  
  cols = brewer.pal(5, name = "Set1")[2:5]
  cols2 = cols[c(1,3)]
  
  plot(smooth.plot(beta.list, Absolute.bias[,1], alpha = alpha1), type = 'l', 
       col = cols[1], lwd = 2, 
       xlab = "",
       ylab = "", 
       xaxt = xaxt,
       ylim = ylim.absolute,
       main = title)
  title(ylab = ylab, line = 2.2)
  title(xlab = xlab, line = 2.2)
  
  lines(smooth.plot(beta.list, abs(coef.save[,8]-coef.save[,1])),
        lwd = 2, lty = 3)
  lines(beta.list, rep(0,length(beta.list)),lwd=2)
  
  
  for(k in 2){
    lines(smooth.plot(beta.list, Absolute.bias[,k], alpha = alpha1),lwd = 2, 
          col = cols[k])
  }
  
  for(k in 3:4){
    lines(smooth.plot(beta.list, Absolute.bias[,k], alpha = alpha2),lwd = 2, 
          col = cols[k])
  }
  
  # estimates from true propensity score 
  for(k in 1:2){
    lines(smooth.plot(beta.list, Absolute.bias[,k+4], alpha = alpha1),lwd = 2, 
          col = cols2[k], lty = 3)
  }
  
  
  
  if(plot.points == TRUE){
    points(beta.list, Absolute.bias[,1],pch = 18, 
           col = cols[1])
    for(k in 2:4){
      points(beta.list, Absolute.bias[,k],pch = 18, 
             col = cols[k])
    }
    for(k in 1:2){
      points(beta.list, Absolute.bias[,k+4],pch = 18, 
             col = cols2[k])
    }
  }
  
  
}

make.plot.notsmooth.coefX = function(coef.save, 
                                     ylim.percent = NULL,
                                     ylim.absolute = NULL,
                                     title = NULL, 
                                     plot.points = TRUE,
                                     xlab = "Coefficient Magnitude",
                                     ylab = "Absolute Bias",
                                     alpha = .8,
                                     beta.list){
  # get percent and absolute bias vs SRS
  SRS = matrix(rep(coef.save[,1],6), nrow = length(beta.list),ncol = 6)
  BS = coef.save[,2:7]
  
  # Percent.bias = abs((BS-SRS)/SRS*100)
  Absolute.bias = abs(BS-SRS)
  
  if(is.null(ylim.absolute)){
    ylim.absolute = c(0,max(Absolute.bias))
  }
  
  cols = brewer.pal(5, name = "Set1")[2:5]
  cols2 = cols[c(1,3)]
  
  plot(beta.list, Absolute.bias[,1], type = 'l', 
       col = cols[1], lwd = 2, 
       xlab = xlab,
       ylab = ylab, ylim = ylim.absolute,
       main = title)
  
  lines(beta.list, abs(coef.save[,8]-coef.save[,1]),
        lwd = 2, lty = 3)
  lines(beta.list, rep(0,length(beta.list)),lwd=2)
  
  
  for(k in 2:4){
    lines(beta.list, Absolute.bias[,k],lwd = 2, 
          col = cols[k])
  }
  
  # estimates from true propensity score 
  for(k in 1:2){
    lines(beta.list, Absolute.bias[,k+4],lwd = 2, 
          col = cols2[k], lty = 3)
  }
  
  
  
  if(plot.points == TRUE){
    points(beta.list, Absolute.bias[,1],pch = 18, 
           col = cols[1])
    for(k in 2:4){
      points(beta.list, Absolute.bias[,k],pch = 18, 
             col = cols[k])
    }
    for(k in 1:2){
      points(beta.list, Absolute.bias[,k+4],pch = 18, 
             col = cols2[k])
    }
  }
  
  
}

make.plot.perc.smooth.coefX = function(coef.save, 
                                       ylim.percent = NULL,
                                       ylim.absolute = NULL,
                                       title = NULL, 
                                       plot.points = TRUE,
                                       xlab = "Coefficient Magnitude",
                                       ylab = "Relative Bias",
                                       alpha = .4,
                                       beta.list){
  # get percent and absolute bias vs SRS
  SRS = matrix(rep(coef.save[,1],6), nrow = length(beta.list),ncol = 6)
  BS = coef.save[,2:7]
  
  # Percent.bias = abs((BS-SRS)/SRS*100)
  Absolute.bias = abs((BS-SRS)/SRS)
  
  if(is.null(ylim.absolute)){
    ylim.absolute = c(0,max(Absolute.bias))
  }
  
  cols = brewer.pal(5, name = "Set1")[2:5]
  cols2 = cols[c(1,3)]
  
  plot(smooth.plot(beta.list, Absolute.bias[,1], alpha = alpha), type = 'l', 
       col = cols[1], lwd = 2, 
       xlab = xlab,
       ylab = ylab, ylim = ylim.absolute,
       main = title)
  
  
  for(k in 2:4){
    lines(smooth.plot(beta.list, Absolute.bias[,k], alpha = alpha),lwd = 2, 
          col = cols[k])
  }
  
  # estimates from true propensity score 
  for(k in 1:2){
    lines(smooth.plot(beta.list, Absolute.bias[,k+4], alpha = alpha),lwd = 2, 
          col = cols2[k], lty = 3)
  }
  
  lines(smooth.plot(beta.list, abs((coef.save[,8]-coef.save[,1])/coef.save[,1])),
        lwd = 2, lty = 3)
  lines(beta.list, rep(0,length(beta.list)),lwd=2)
  
  if(plot.points == TRUE){
    points(beta.list, Absolute.bias[,1],pch = 18, 
           col = cols[1])
    for(k in 2:4){
      points(beta.list, Absolute.bias[,k],pch = 18, 
             col = cols[k])
    }
    # for(k in 1:2){
    #   points(beta.list, Absolute.bias[,k+4],pch = 18, 
    #          col = cols2[k])
    # }
  }
  
  
}

colnames = c("SRS", paste("Biased:",c("No weighting", "Weighted prop model", "Weighted outcome model", "Both weighted",
                                      "No weighting - True prop", "Weighted outcome - True prop")),
             "SRS weighted")

cols = brewer.pal(5, name = "Set1")[2:5]
cols2 = cols[c(1,3)]

path = "~/Dropbox/UCI/Research/Prediction Assessment/weighted_propensity_scores/simstudy/"

#### Vary correlation ####
# results from simstudy_vary_corr.R
load(paste0(path,"coefvary_corr.RData"))
coefvary_mvn = do.call("rbind",simresults)


pdf(paste0(path,"figures/varycorr.pdf"),
    width = 6.5, height = 5)

par(mfrow = c(1,1), 
    mar = c(3.1, 3.1, 1.6, .6), 
    oma = c(.3,.3,.5,6.3), xpd = NA)

c.list = c(.05*(1:19),.96,.97,.98,.99,.999,1)

make.plot.perc.smooth.coefX(coefvary_mvn,#[27:51,],
                            #ylim.absolute = c(0,2),
                            plot.points = FALSE,
                            title = "Vary correlation",
                            xlab = "", 
                            ylab = "",
                            beta.list = c.list,
                            alpha = .6)

title(xlab = expression(paste("Correlation: ",rho)), line = 2.2)
title(ylab = "Relative Bias (vs SRS)", line = 2.2)

legend(1.04,.55,
       legend = c("SRS\nEst. prop. score","SRS\nTrue prop. score",
                  paste0("Conv. Sample: ", c("\nNo weights\nEst prop. score",
                                      "\nNo weights\nTrue prop. score",
                                      "\nWeight prop. model\nEst. prop. score", 
                                      "\nWeight outcome model\nEst. prop. score",
                                      "\nWeight outcome model\nTrue prop. score", 
                                      "\nWeight both models\nEst. prop score"))), 
       col = c("black","black",cols[1], cols[1],cols[2:3],cols[3:4],cols[4:5]), lwd = rep(2,8),
       lty = c(1,3,1,3,1,1,3,1),
       bty = 'n', cex = .6,
       y.intersp=2.1)

dev.off()

#### Vary coefficients in outcome model, correlation = 1 ####
# code from ~/obernste/propscore/simstudy/simstudy_vary_TK.R,
# /simstudy_vary_T.R, /simstudy_vary_X2.R, /simstudy_vary_X3.R
# /simstudy_vary_X4.R, /simstudy_vary_X5.R

pdf(paste0(path,"figures/varycoef_cor1.pdf"),
    width = 6.5, height = 5)

par(mfrow = c(2,3), 
    mar = c(3.1, 3.1, 1.6, .6), 
    oma = c(.3,.3,2.5,10), xpd = NA)

beta.list = (-10:10)/2

# results from simstudy_vary_T.R
load(paste0(path,"coefvary_T_1.RData"))
coefvary_T = do.call("rbind",simresults)
make.plot.smooth.coefX(coefvary_T,
                       ylim.absolute = c(0,3),
                       plot.points = FALSE,
                       title = expression(paste("Vary ",delta[A])),
                       ylab = "Absolute Bias (vs SRS)",
                       beta.list = beta.list,
                       xlab = expression(delta[A]))

# results from simstudy_vary_TK.R
load(paste0(path,"coefvary_TK_1.RData"))
coefvary_TK = do.call("rbind",simresults)
make.plot.smooth.coefX(coefvary_TK,
                       ylim.absolute = c(0,3),
                       plot.points = FALSE,
                       alpha1 = .4,
                       alpha2 = .9,
                       title = expression(paste("Vary ",delta[AK])),
                       ylab = "Absolute Bias (vs SRS)",
                       beta.list = beta.list,
                       xlab = expression(delta[AK]))

# results from simstudy_vary_X2.R
load(paste0(path,"coefvary_X2_1.RData"))
coefvary_X2 = do.call("rbind",simresults)
make.plot.smooth.coefX(coefvary_X2,
                       ylim.absolute = c(0,3),
                       plot.points = FALSE,
                       title = expression(paste("Vary ",delta[2])),
                       ylab = "Absolute Bias (vs SRS)",
                       beta.list = beta.list,
                       xlab = expression(delta[2]))


legend(5.2,3.5,
       legend = c("SRS\nEst. prop. score","SRS\nTrue prop. score",
                  paste0("Conv. Sample: ", c("\nNo weights\nEst prop. score",
                                             "\nNo weights\nTrue prop. score",
                                             "\nWeight prop. model\nEst. prop. score",
                                             "\nWeight outcome model\nEst. prop. score",
                                             "\nWeight outcome model\nTrue prop. score",
                                             "\nWeight both models\nEst. prop score"))),
       col = c("black","black",cols[1], cols[1],cols[2:3],cols[3:4],cols[4:5]), lwd = rep(2,8),
       lty = c(1,3,1,3,1,1,3,1),
       bty = 'n', cex = .94,
       y.intersp=2.1)

# results from simstudy_vary_X3.R
load(paste0(path,"coefvary_X3_1.RData"))
coefvary_X3 = do.call("rbind",simresults)
make.plot.smooth.coefX(coefvary_X3,
                       ylim.absolute = c(0,3),
                       plot.points = FALSE,
                       title = expression(paste("Vary ",delta[3])),
                       ylab = "Absolute Bias (vs SRS)",
                       beta.list = beta.list,
                       xlab = expression(delta[3]))

# results from simstudy_vary_T.R
load(paste0(path,"coefvary_X4_1.RData"))
coefvary_X4 = do.call("rbind",simresults)
make.plot.smooth.coefX(coefvary_X4,
                       ylim.absolute = c(0,3),
                       plot.points = FALSE,
                       title = expression(paste("Vary ",delta[4])),
                       ylab = "Absolute Bias (vs SRS)",
                       beta.list = beta.list,
                       xlab = expression(delta[4]))

# results from simstudy_vary_X5.R
load(paste0(path,"coefvary_X5_1.RData"))
coefvary_X5 = do.call("rbind",simresults)
make.plot.smooth.coefX(coefvary_X5,
                       ylim.absolute = c(0,3),
                       plot.points = FALSE,
                       title = expression(paste("Vary ",delta[5])),
                       ylab = "Absolute Bias (vs SRS)",
                       beta.list = beta.list,
                       xlab =  expression(delta[5]))


title("Vary Coefficients in the True Outcome Model:\n K is Observed", 
      outer = TRUE, line = 0.4)

dev.off()

#### comparing varying delta_TK with different correlations
beta.list = (-10:10)/2

corr = .95
Sigma = matrix(c(1,corr,corr,1),nrow = 2)
cont = MASS::mvrnorm(1000000, mu = c(0,0), Sigma = Sigma)
K.cont = cont[,1]
x1.cont = cont[,2]

K =  ifelse(K.cont > 0,1,0)
x1 = ifelse(x1.cont > 0,1,0)
misclass_95 = rnd(mean(K!=x1),2)

corr = .9
Sigma = matrix(c(1,corr,corr,1),nrow = 2)
cont = MASS::mvrnorm(1000000, mu = c(0,0), Sigma = Sigma)
K.cont = cont[,1]
x1.cont = cont[,2]

K =  ifelse(K.cont > 0,1,0)
x1 = ifelse(x1.cont > 0,1,0)
misclass_9 = rnd(mean(K!=x1),2)

pdf(paste0(path,"figures/vary_TK_and_corr.pdf"),
    width = 6.5, height = 2.5)

par(mfrow = c(1,3), 
    mar = c(3.1, 3.1, 1.6, .6), 
    oma = c(.3,.3,2.5,7.3), xpd = NA)

# results from simstudy_vary_TK.R
load(paste0(path,"coefvary_TK_1.RData"))
coefvary_TK = do.call("rbind",simresults)
make.plot.smooth.coefX(coefvary_TK,
                       ylim.absolute = c(0,2),
                       plot.points = FALSE,
                       alpha1 = .4,
                       alpha2 = .6,
                       title = expression(paste(rho, " = 1, Misclass: 0")),
                       ylab = "Absolute Bias (vs SRS)",
                       beta.list = beta.list,
                       xlab = expression(delta[AK]))

# results from simstudy_vary_TK.R but need to change corr = .95
load(paste0(path,"coefvary_TK_95.RData"))
coefvary_TK_95 = do.call("rbind",simresults)
make.plot.smooth.coefX(coefvary_TK_95,
                       ylim.absolute = c(0,2),
                       plot.points = FALSE,
                       alpha1 = .4,
                       alpha2 = .6,
                       title = expression(paste(rho, " = .95, Misclass: .10")),
                       ylab = "Absolute Bias (vs SRS)",
                       beta.list = beta.list,
                       xlab = expression(delta[AK]))

# results from simstudy_vary_TK.R but need to change corr = .9
load(paste0(path,"coefvary_TK_9.RData"))
coefvary_TK_9 = do.call("rbind",simresults)
make.plot.smooth.coefX(coefvary_TK_9,
                       ylim.absolute = c(0,2),
                       plot.points = FALSE,
                       alpha1 = .4,
                       alpha2 = .6,
                       title = expression(paste(rho, " = .90, Misclass: .14",)),
                       ylab = "Absolute Bias (vs SRS)",
                       beta.list = beta.list,
                       xlab = expression(delta[AK]))

legend(6,2.9,
       legend = c("SRS\nEst. prop. score","SRS\nTrue prop. score",
                  paste0("Conv. Sample: ", c("\nNo weights\nEst prop. score",
                                             "\nNo weights\nTrue prop. score",
                                             "\nWeight prop. model\nEst. prop. score",
                                             "\nWeight outcome model\nEst. prop. score",
                                             "\nWeight outcome model\nTrue prop. score",
                                             "\nWeight both models\nEst. prop score"))),
       col = c("black","black",cols[1], cols[1],cols[2:3],cols[3:4],cols[4:5]), lwd = rep(2,8),
       lty = c(1,3,1,3,1,1,3,1),
       bty = 'n', cex = .6,
       y.intersp=2.1)

title(expression(paste("Vary ",delta[AK], " and ",rho)), outer = TRUE, line = 1, text = 2)

dev.off()


#### varying coefficients in the propensity score model ####
library(patchwork)
library(latticeExtra)
library(viridis)

beta.list = .1*(5:20)

load(paste0(path,"coefvary_propmodel_2.RData"))
prop_2 = simresults

load(paste0(path,"coefvary_propmodel_3.RData"))
prop_3 = simresults


data1 <- expand.grid(K0=beta.list, K1=beta.list)
data1$bias_2<- unlist(lapply(prop_2,function(x){abs(x[4]-x[1])}))
data1$bias_3<- unlist(lapply(prop_3,function(x){abs(x[4]-x[1])}))
data1$model = "Weighted Outcome Model"

# weighting both models
data2 <- expand.grid(K0=beta.list, K1=beta.list)
data2$bias_2<- unlist(lapply(prop_2,function(x){abs(x[5]-x[1])}))
data2$bias_3<- unlist(lapply(prop_3,function(x){abs(x[5]-x[1])}))
data2$model = "Weighted Propensity and Outcome Models"

data.comb = rbind(data1, data2)
data.comb$model = as.factor(data.comb$model)

jpeg(paste0(path,"figures/vary_promodel_X2.jpeg"),
     width = 8.5, height = 5, 
     units = "in", res = 1000)

levelplot(bias_2 ~ K0 * K1 | model, data.comb, 
          panel = panel.2dsmoother, n = 200,
          xlab = expression(paste(alpha["02"],": Exp. Coef. on (1-K)*",X[2], " in Propensity Model")),
          ylab = expression(paste(alpha["12"],": Exp. Coef. on K*",X[2], " in Propensity Model")),
          main = "Absolute Bias in Treatment Effect Estimate Relative to SRS",
          col.regions = viridis(100), 
          contour = TRUE,
          at = .025*(-2:26), # functions as a zlim
          scales = list(alternating = 1, tck = c(1,0)), 
          strip = strip.custom(bg = "grey90"))

dev.off()

jpeg(paste0(path,"figures/vary_promodel_X3.jpeg"),
     width = 8.5, height = 5, 
     units = "in", res = 1000)

levelplot(bias_3 ~ K0 * K1 | model, data.comb, 
          panel = panel.2dsmoother, n = 200,
          xlab = expression(paste(alpha["03"],": Exp. Coef. on (1-K)*",X[3], " in Propensity Model")),
          ylab = expression(paste(alpha["13"],": Exp. Coef. on K*",X[3], " in Propensity Model")),
          main = "Absolute Bias in Treatment Effect Estimate Relative to SRS",
          col.regions = viridis(100), 
          contour = TRUE, 
          at = .025*(-2:26), # functions as a zlim
          scales = list(alternating = 1, tck = c(1,0)), 
          strip = strip.custom(bg = "grey90"))

dev.off()

