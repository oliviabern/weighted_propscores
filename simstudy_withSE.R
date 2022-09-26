# Update simulation scenario 4/13/22

# fix coefficients and compare standard error estimates

library("splines")
library("tidyverse")
library("survey")
library("snowfall")
library("rlecuyer")
library("estweight")

n.pop = 50000
n = 2000
n.sim = 500

# source("~/Documents/UCI/18W/211/Stat211Functions.R")
# source("~/Dropbox/helper.R")

## function from Dan Gillen's 211 code for estimating the robust sandwich se estimate
#####	robust.se.lm() is a function to compute the Huber-White sandwich variance estimator
#####	for the linear regression model
#####	
##
robust.se.lm <- function( model) { 
  s <- summary( model) 
  X <- model.matrix( model )
  sandwich.cov <- robust.vcov.lm( model )
  sand.se <- sqrt( diag( sandwich.cov )) 
  t <- model$coefficients/sand.se
  p <- 2*pt( -abs( t ), dim(X)[1]-dim(X)[2] ) 
  ci95.lo <- model$coefficients - qt( .975, dim(X)[1]-dim(X)[2] ) * sand.se
  ci95.hi <- model$coefficients + qt( .975, dim(X)[1]-dim(X)[2] ) * sand.se
  rslt <- cbind( model$coefficients, sand.se, ci95.lo, ci95.hi, t, p ) 
  dimnames(rslt)[[2]] <- c( dimnames( s$coefficients )[[2]][1], "Robust SE", "ci95.lo", "ci95.hi", dimnames( s$coefficients )[[2]][3:4] ) 
  rslt 
} 

##
##### 	Compute robust (sandwich) variance-covariance estimate for a LM
##
robust.vcov.lm <- function( lm.obj ){
  X <- model.matrix( lm.obj )
  eps <- lm.obj$residuals
  robust.cov <- solve( t(X)%*%X ) %*%( t(X) %*% diag(eps^2) %*% X ) %*% solve( t(X)%*%X )
  dimnames( robust.cov ) <- dimnames( vcov(lm.obj) )
  return( robust.cov )
}

# round variable and set the numebr of trailing 0's
rnd = function(x, n.dig){format(round(x, digits = n.dig), nsmall = n.dig)}

expit <- function(x){exp(x)/(1+exp(x))}
logit = function(x){log(x/(1-x))}

outcome_family = stats::gaussian

p.b = 0

coef.srs <- coef.nowt <- coef.wt.prop <- coef.wt.outcome <-coef.wt2 <- coef.wt.outcome.trueprop <- coef.nowt.trueprop <- coef.srs.trueprop <- matrix(NA, nrow = n.sim, ncol = 3)
se.srs <- se.nowt <- se.wt.prop  <- se.wt.outcome <- se.wt2.prop  <- se.wt2.db <- se.wt2.naive <- matrix(NA, nrow = n.sim, ncol = 3)

beta.t = 1
beta.x2 = 1.5
beta.x3 = -2
beta.x4 = -1
beta.x5 = 1.5
beta.t.k = 3
beta.1_k.x2 = 1.3
beta.k.x2.2 = 2
beta.1_k.x3 = .4#~.83
beta.k.x3.2 = 1.5

corr = 1

set.seed(1364401596)

start = Sys.time()

for(i in 1:n.sim){
  
  if(round(i/50)*50 == i){print(i)}
  
  # generate biased sample
  Sigma = matrix(c(1,corr,corr,1),nrow = 2)
  
  cont = MASS::mvrnorm(n.pop, mu = c(p.b, p.b), Sigma = Sigma)
  K.cont = cont[,1]
  x1.cont = cont[,2]
  
  K =  ifelse(K.cont > 0,1,0)
  x1 = ifelse(x1.cont > 0,1,0)
  
  
  # prob.sample = pnorm(K)
  prob.sample = .6*K + .2
  
  # get biased sample
  BSamp = sample(1:n.pop, size = n, prob = prob.sample)
  x1.b = x1[BSamp]
  K.b = K[BSamp]
  
  # get SRS
  SRS = sample(1:n.pop, size = n)
  x1.s = x1[SRS]
  K.s = K[SRS]
  
  # get remaining covariates for biased sample
  x2.b = rnorm(n,0,2)
  x3.b = rnorm(n,0,2)
  pi.b = expit((log(beta.1_k.x2)*x2.b 
                + log(beta.1_k.x3)*x3.b)*(1-K.b) + 
                 (log(beta.k.x2.2)*x2.b 
                  + log(beta.k.x3.2)*x3.b)*(K.b))
  T.b = rbinom(n,1,pi.b) 
  x4.b = rnorm(n,0,1)
  x5.b = rnorm(n,0,1)
  y.b = rnorm(n, mean = (0 + beta.t*T.b + beta.t.k*T.b*K.b +
                           beta.x2*x2.b + beta.x3*x3.b^2 + 
                           beta.x4*x4.b + beta.x5*x5.b^3), sd = 1) 
  
  # generate SRS
  x2.s = rnorm(n,0,2)
  x3.s = rnorm(n,0,2)
  pi.s = expit((log(beta.1_k.x2)*x2.s 
                + log(beta.1_k.x3)*x3.s)*(1-K.s) + 
                 (log(beta.k.x2.2)*x2.s 
                  + log(beta.k.x3.2)*x3.s)*(K.s))
  T.s = rbinom(n,1,pi.s) 
  x4.s = rnorm(n,0,1)
  x5.s = rnorm(n,0,1)
  y.s = rnorm(n, mean = (0 + beta.t*T.s + beta.t.k*T.s*K.s +
                           beta.x2*x2.s + beta.x3*x3.s^2 + 
                           beta.x4*x4.s + beta.x5*x5.s^3), sd = 1) 
  
  # # estimate sampling weights for biased sample
  biased = c(rep(1,n),rep(0,n))
  data.comb = data.frame(
    biased = biased,
    K.comb = c(K.b, K.s),
    x1.comb = c(x1.b, x1.s),
    x2.comb = c(x2.b, x2.s),
    x3.comb = c(x3.b, x3.s),
    x4.comb = c(x4.b, x4.s),
    x5.comb = c(x5.b, x5.s)
  )
  
  # set scope for stepwise AIC
  vars = colnames(data.comb)[c(3:7)]
  
  form_min = stats::as.formula("biased ~ 1")
  form_max = stats::as.formula(paste0("biased ~",paste0(c(
    vars, paste0("I(",vars,"^2)")
  ),collapse = "+")))
  fit_min = stats::glm(formula = form_min, family = stats::binomial, data = data.comb)
  forward = stats::step(fit_min,scope=list(lower=form_min,upper=form_max),
                        direction="forward", trace = 0)
  
  # fit selected model
  estwt_form = stats::formula(forward)
  
  estwt.fit = stats::glm(estwt_form,
                         family = binomial, data = data.comb)
  
  prob.bias = fitted(estwt.fit, "response")       # est. probability biased == 1
  htweight_unnorm = (1-prob.bias)/prob.bias       # HT weight
  htweight = (htweight_unnorm/sum(htweight_unnorm))[1:n]
  
  
  #### Weighted propensity score model #####
  
  # estimate propensity score (with weights)
  data.use = data.frame(T.b, K.b, 
                        x1.b, x2.b, 
                        x3.b, x4.b, x5.b,
                        y.b, htweight)
  
  data.use.s = data.frame(T.s, K.s,
                          x1.s, x2.s,
                          x3.s, x4.s, x5.s,
                          y.s) 
  
  ## set scope for weighted stepwise AIC
  vars_prop = colnames(data.use)[c(3:7)]
  
  form_min_prop = stats::as.formula("T.b ~ 1")
  form_max_prop = stats::as.formula(paste0("T.b ~",paste0(c(
    vars_prop, 
    paste0(vars_prop[1],":",vars_prop),
    paste0(vars_prop[2],":",vars_prop),
    paste0(vars_prop[3],":",vars_prop),
    paste0(vars_prop[4],":",vars_prop),
    paste0(vars_prop[5],":",vars_prop)
  ),collapse = "+")))
  
  fit_min = do.call("svyglm", args = list(form_min_prop,
                                          design = survey::svydesign(ids = ~0, weights = data.use$htweight, data = data.use),
                                          family = stats::quasibinomial))
  
  forward = stats::step(fit_min,scope=list(lower=form_min_prop,upper=form_max_prop),
                        direction="forward", trace = 0)
  # fit selected model (weighted)
  form.b.weighted = stats::formula(forward)
  
  
  ## select model for unweighted stepwise AIC
  fit_min_unwt = glm(form_min_prop,
                     data = data.use,
                     family = stats::binomial)
  forward_unwt = stats::step(fit_min_unwt,scope=list(lower=form_min_prop,upper=form_max_prop),
                             direction="forward", trace = 0)
  # fit selected model (unweighted)
  form.b.unwt = stats::formula(forward_unwt)
  
  # ## set scope for stepwise AIC in SRS
  vars_prop_s = colnames(data.use.s)[c(3:7)]
  
  form_min_prop_s = stats::as.formula("T.s ~ 1")
  form_max_prop_s = stats::as.formula(paste0("T.b ~",paste0(c(
    vars_prop_s, 
    paste0(vars_prop_s[1],":",vars_prop_s),
    paste0(vars_prop_s[2],":",vars_prop_s),
    paste0(vars_prop_s[3],":",vars_prop_s),
    paste0(vars_prop_s[4],":",vars_prop_s),
    paste0(vars_prop_s[5],":",vars_prop_s)
  ),collapse = "+")))
  fit_min_s = glm(form_min_prop_s, data = data.use.s,
                  family = stats::binomial)
  forward_s = stats::step(fit_min_s,scope=list(lower=form_min_prop_s,upper=form_max_prop_s),
                          direction="forward", trace = 0)
  
  # fit selected model (weighted)
  form.s = stats::formula(forward_s)
  estprop.nowt = glm(form.b.unwt,data = data.use, family = binomial)
  estprop = svyglm(form.b.weighted, 
         design = survey::svydesign(ids = ~0, weights = data.use$htweight, data = data.use),
         family = stats::quasibinomial)
  estprop.srs = glm(form.s,data = data.use.s, family = binomial)
  
  
  # try mean centering (removed mean centering)
  p = fitted(estprop,"response")
  data.use$propscore = (p - mean(p))/sd(p)
  p2 = fitted(estprop.nowt, "response")
  data.use$propscore.nowt = (p2 - mean(p2))/sd(p2)
  p3 = fitted(estprop.srs)
  data.use.s$propscore.srs  = (p3 - mean(p3))/sd(p3)
  
  data.use$T.b = data.use$T.b - mean(data.use$T.b)
  data.use.s$T.s = data.use.s$T.s - mean(data.use.s$T.s)
  
  
  # fit outcome model (weights X2)
  # data.use$T.b.meancent = data.use$T.b - mean(data.use$T.b);
  outcome_wts2 = svyglm(y.b ~ T.b + propscore, 
                        design = survey::svydesign(ids = ~0, weights = data.use$htweight, data = data.use),
                        family = stats::gaussian)
  coef.wt2[i,] = outcome_wts2$coefficients
  
  # design based estimate
  se.wt2.prop[i,] = treatmentSE(estwt_fit = estwt.fit, 
                           estprop_fit = estprop, 
                           fit_outcome = outcome_wts2, 
                           biased = biased,
                           outcome_family = outcome_family)
  
  se.wt2.db[i,] = summary(outcome_wts2)$coef[,2]
  
  # naive (non sandwich) estimate
  outcome_naive = glm(y.b ~ T.b + propscore, 
                      weights = data.use$htweight,
                      data = data.use,
                      family = stats::gaussian)
  se.wt2.naive[i,] = summary(outcome_naive)$coef[,2]
  
  # fit outcome model (weights only in prop model)
  outcome = glm(y.b ~ T.b + propscore, 
                data = data.use,
                family = stats::gaussian)
  coef.wt.prop[i,] = outcome$coefficients
  
  se.wt.prop[i,] = robust.se.lm(outcome)[,2]
  
  # fit outcome model (weights only in outcome model)
  outcome = svyglm(y.b ~ T.b + propscore.nowt, 
                   design = survey::svydesign(ids = ~0, weights = data.use$htweight, data = data.use),
                   family = stats::gaussian)
  coef.wt.outcome[i,] = outcome$coefficients
  
  se.wt.outcome[i,] = summary(outcome_naive)$coef[,2]
  
  # fit outcome model (no weighting)
  outcome2 = glm(y.b ~ T.b + propscore.nowt, 
                 data = data.use,
                 family = stats::gaussian)
  coef.nowt[i,] = outcome2$coefficients
  
  se.nowt[i,] = robust.se.lm(outcome2)[,2]
  
  # fit outcome model in SRS
  # data.use.s$T.s.meancent = data.use.s$T.s - mean(data.use.s$T.s);
  outcome = glm(y.s ~ T.s + propscore.srs, 
                data = data.use.s,
                family = stats::gaussian)
  coef.srs[i,] = outcome$coefficients
  
  se.srs[i,] = robust.se.lm(outcome)[,2]

  
}

end = Sys.time()
  
end-start

coef.out = rbind(apply(coef.srs,2,mean),
                 apply(coef.nowt, 2, mean),
                 apply(coef.wt.prop, 2, mean),
                 apply(coef.wt.outcome, 2, mean),
                 apply(coef.wt2, 2, mean)#,
                 # apply(coef.nowt.trueprop, 2, mean),
                 # apply(coef.wt.outcome.trueprop, 2, mean),
                 # apply(coef.srs.trueprop, 2, mean)
                 )
colnames(coef.out) = c("Intercept","Treatment","Prop Score")
rownames(coef.out) = c("SRS",paste("Biased:",
                                   c("No weights", "Wt prop model", 
                                     "Wt outcome model","Wt both models")))

se.analytic = rbind(apply(se.srs,2,mean),
                    apply(se.nowt, 2, mean),
                    apply(se.wt.prop, 2, mean),
                    apply(se.wt.outcome, 2, mean),
                    apply(se.wt2.prop, 2, mean))
dimnames(se.analytic) = dimnames(coef.out)

se.empirical = rbind(apply(coef.srs,2,sd),
                     apply(coef.nowt, 2, sd),
                     apply(coef.wt.prop, 2, sd),
                     apply(coef.wt.outcome, 2, sd),
                     apply(coef.wt2, 2, sd))

ci.lo.emp = rbind(apply(coef.srs,2,function(x){quantile(x,.025)}),
                  apply(coef.nowt, 2, function(x){quantile(x,.025)}),
                  apply(coef.wt.prop, 2, function(x){quantile(x,.025)}),
                  apply(coef.wt.outcome, 2, function(x){quantile(x,.025)}),
                  apply(coef.wt2, 2, function(x){quantile(x,.025)}))

ci.hi.emp = rbind(apply(coef.srs,2,function(x){quantile(x,.975)}),
                  apply(coef.nowt, 2, function(x){quantile(x,.975)}),
                  apply(coef.wt.prop, 2, function(x){quantile(x,.975)}),
                  apply(coef.wt.outcome, 2, function(x){quantile(x,.975)}),
                  apply(coef.wt2, 2, function(x){quantile(x,.975)}))

dimnames(se.empirical) = dimnames(coef.out)

# compare the proposed, design based, and naive estimators
# in the doubly weighted method
se.comparemethods = rbind(apply(coef.wt2, 2, sd),
                          apply(se.wt2.naive, 2, mean),
                          apply(se.wt2.db, 2, mean),
                          apply(se.wt2.prop, 2, mean))
colnames(se.comparemethods) = c("Intercept","Treatment","Prop Score")
rownames(se.comparemethods) = c("Empirical", "Naive (non-sandwich)", "Design Based", "Proposed")
round(se.comparemethods,2)

summary(coef.wt2)
hist(coef.wt2[,3])

# I think there must be some extreme values of coefficients
# probably coming from extreme weights or small subpopulations

# #### forest plot of point estimates and analytic variance estimates ####
# se.bars.vert = function(ci.lo, ci.hi, y, col, lwd = 3){
#   arrows(ci.lo, y, ci.hi, y, length=0.05, angle=90, code=3, col = col,lwd = lwd)
# }
# 
# # Generate point estimates and standard errors
# pt.est = coef.out[,2]
# se = 1.96*se.analytic[,2]
# 
# # calculate CI
# ci.lo = pt.est - se
# ci.hi = pt.est + se
# 
# # create string for table
# string = paste0(round(pt.est,2),
#                 " (",
#                 round(ci.lo,2),
#                 ", ", 
#                 round(ci.hi,2),")")
# 
# # Pick spots on y-axis for plotting
# y = length(pt.est):1
# 
# # change plotting parameters to leave room for text
# # on the left side 
# par(oma = c(0,13,0,0) + 0.1, # adjust the second entry in the oma command to shift left side margin
#     mar = c(3,0,4,1) + 0.1, ps = 10)
# 
# 
# plot(pt.est,y, 
#      xlim = range(c(ci.lo,ci.hi)),
#      ylim = c(min(y)-1,max(y)+1), 
#      yaxt = 'n', cex = 1.5,
#      pch = 18, ylab = "", xlab = "",
#      main = "Estimated Average Treatment Effect")
# se.bars.vert(ci.lo, ci.hi, y = y, col = 'black')
# abline(v = pt.est[1], lty = 3, lwd = 2, col = 'grey')
# 
# # add text 
# mtext(c("Mean (95% CI)",string), 2, line = 1, # adjust the line value to shift text left and right
#       las = 2, at = c(max(y)+1,y), adj = 1)
# models = c("SRS",paste("Biased:\n",
#                        c("No weights", "Wt prop model", 
#                          "Wt outcome model","Wt both models")))
# mtext(models, 
#       2, line = 7, font = 2,
#       las = 2, at = y, adj = 1)
# 
# #### forest plot of point estimates and empirical variance estimates
# pt.est = coef.out[,2]
# # se = 1.96*se.analytic[,2]
# 
# # calculate CI
# ci.lo = ci.lo.emp[,2]
# ci.hi = ci.hi.emp[,2]
# 
# # create string for table
# string = paste0(round(pt.est,2),
#                 " (",
#                 round(ci.lo,2),
#                 ", ", 
#                 round(ci.hi,2),")")
# 
# plot(pt.est,y, 
#      xlim = range(c(ci.lo,ci.hi)),
#      ylim = c(min(y)-1,max(y)+1), 
#      yaxt = 'n', cex = 1.5,
#      pch = 18, ylab = "", xlab = "",
#      main = "Estimated Average Treatment Effect")
# se.bars.vert(ci.lo, ci.hi, y = y, col = 'black')
# abline(v = pt.est[1], lty = 3, lwd = 2, col = 'grey')
# 
# # add text 
# mtext(c("Mean (95% CI)",string), 2, line = 1, # adjust the line value to shift text left and right
#       las = 2, at = c(max(y)+1,y), adj = 1)
# models = c("SRS",paste("Biased:\n",
#                        c("No weights", "Wt prop model", 
#                          "Wt outcome model","Wt both models")))
# mtext(models, 
#       2, line = 7, font = 2,
#       las = 2, at = y, adj = 1)
# 

### forest plot of point estimates +/- 2*empirical SD plus analytic estimates for doubly weighted model ####
se.bars.vert = function(ci.lo, ci.hi, y, col, lwd = 3){
  arrows(ci.lo, y, ci.hi, y, length=0.05, angle=90, code=3, col = col,lwd = lwd)
}

# Generate point estimates and standard errors
pt.est = coef.out[,2]
se = 1.96*se.empirical[,2]

# calculate CI
ci.lo = pt.est - se
ci.hi = pt.est + se

# get CI's for doubly weighted mode

# create string for table
string = paste0(rnd(pt.est,2),
                " (",
                rnd(ci.lo,2),
                ", ",
                rnd(ci.hi,2),")")

# Pick spots on y-axis for plotting
y = length(pt.est):1

pdf("~/Dropbox/UCI/Research/Prediction Assessment/weighted_propensity_scores/simstudy/figures/simstudy_withSE_cor1.pdf",
    width = 6.5, height = 3.5)

# change plotting parameters to leave room for text
# on the left side
par(oma = c(0,13,0,9.3) + 0.1, # adjust the second entry in the oma command to shift left side margin
    mar = c(3,0,1,1) + 0.1, ps = 10, xpd = NA)

plot(pt.est,y,
     #xlim = range(c(ci.lo,ci.hi)),
     xlim = c(0,5.8),
     ylim = c(min(y)-1,max(y)+1),
     yaxt = 'n', cex = 1.5, 
     pch = 18, ylab = "", xlab = "Mean (95% CI)",
     main = "")
lines(rep(pt.est[1],2),c(min(y)-1.2,max(y)+1.2), lty = 3, lwd = 2, col = 'grey')
se.bars.vert(ci.lo, ci.hi, y = y, col = 'black')


# add text
mtext("Mean Estimate (95% CI)", 1,
      at = 2.5, line = 2)

mtext(c("      Average Estimated\n Treatment Effect (95% CI)"), 2, line = 11, # adjust the line value to shift text left and right
      las = 2, at = max(y)+1.2, adj = 0, font = 2, cex = 1.1)
mtext(string, 2, line = 0.5, 
      las = 2, at = y, adj = 1, font = 1)

mtext(c("Simple Random Sample","Convenience Sample"),
      2, line = 13, font = 2,
      las = 2, at = c(5.5, 4.5), adj = 0)

models = c("No weights",
           "Weight prop. model",
           "Weight outcome model",
           "Weight both models")
mtext(models,
      2, line = 13, font = 1,
      las = 2, at = 4:1, adj = 0)

# add table of SE estimates on the right

mtext(c("Empirical","Proposed","Design Based", "Naive"),
      4, line = 1.5, font = 1, cex = 1.1,
      las = 2, at = (8:5)/2, adj = 0)

mtext(rnd(se.comparemethods[c(1,4:2),2],2),
      4, line = 8.25, font = 1, cex = 1.1,
      las = 2, at = (8:5)/2, adj = 0)

left.x = 6.3 # was 5.5
mid.x = 11.2 # was 9.6
right.x = 13.4 # was 11.6

segments(x0 = left.x, y0 = 5.5, y1 = 2)
segments(x0 = mid.x, y0 = 4.5, y1 = 2)
segments(x0 = right.x, y0 = 5.5, y1 = 2)

segments(x0 = left.x, x1 = right.x, y0 = 5.5)
segments(x0 = left.x, x1 = right.x, y0 = 2)
segments(x0 = left.x, x1 = right.x, y0 = 4.5)

mtext("SE Estimate for Estimate\nWith Both Models Weighted",
      4, line = .8, cex = 1.1, 
      las = 2, at = 5, adj = 0)

dev.off()

# add table to the right
# match left side of table with applied example

