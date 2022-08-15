# Varying correlation

# adding in estimates with K in the scope

# Specify scenario
scenario = "vary_X3_1"

library("splines")
library("tidyverse")
library("survey")
library("snowfall")
library("rlecuyer")
library("rje")

n = 2000
n.sim = 1000

corr = 1

# path = "~/Dropbox/UCI/Research/Prediction Assessment/weighted_propensity_scores/simstudy/"
path = "/home/obernste/propscore/simstudy/"

expit <- function(x){exp(x)/(1+exp(x))}
logit = function(x){log(x/(1-x))}

outcome_family = stats::gaussian

p.b = 0 

coef.srs <- coef.nowt <- coef.wt.prop <- coef.wt.outcome <-coef.wt2 <- coef.wt.outcome.trueprop <- coef.nowt.trueprop <- coef.srs.trueprop <- coef.wt2.Kprop <- coef.wt.outcome.Kprop <- coef.wt.outcome.Kwt <- coef.wt2.Kwt <- matrix(NA, nrow = n.sim, ncol = 3)

beta.t = 1
beta.x2 = 1.5
# beta.x3 = -2
beta.x4 = -1
beta.x5 = 1.5
beta.t.k = 3
beta.1_k.x2 = 1.3
beta.k.x2 = 2
beta.1_k.x3 = .4
beta.k.x3 = 1.5


beta.list = (-10:10)/2

n.pop = 10000

sim = function(beta){
  
  beta.x3 = beta 
  
  for(i in 1:n.sim){
    
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
                   (log(beta.k.x2)*x2.b 
                    + log(beta.k.x3)*x3.b)*(K.b))
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
                   (log(beta.k.x2)*x2.s 
                    + log(beta.k.x3)*x3.s)*(K.s))
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
    
    # # set scope for stepwise AIC
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
    estprop = glm(form.b.weighted,data = data.use, weights = htweight, family = quasibinomial)
    estprop.srs = glm(form.s,data = data.use.s, family = binomial)
    
    
    data.use$propscore = fitted(estprop,"response")
    data.use$propscore.nowt = fitted(estprop.nowt, "response")
    data.use.s$propscore.srs  = fitted(estprop.srs)
    
    
    # fit outcome model (weights X2)
    outcome = svyglm(y.b ~ T.b + propscore, 
                     design = survey::svydesign(ids = ~0, weights = data.use$htweight, data = data.use),
                     family = stats::gaussian)
    coef.wt2[i,] = outcome$coefficients
    
    # fit outcome model (weights only in prop model)
    outcome = glm(y.b ~ T.b + propscore, 
                  data = data.use,
                  family = stats::gaussian)
    coef.wt.prop[i,] = outcome$coefficients
    
    
    # fit outcome model (weights only in outcome model)
    outcome = svyglm(y.b ~ T.b + propscore.nowt, 
                     design = survey::svydesign(ids = ~0, weights = data.use$htweight, data = data.use),
                     family = stats::gaussian)
    coef.wt.outcome[i,] = outcome$coefficients
    
    # fit outcome model (no weighting)
    outcome2 = glm(y.b ~ T.b + propscore.nowt, 
                   data = data.use,
                   family = stats::gaussian)
    coef.nowt[i,] = outcome2$coefficients
    
    # estimate propensity score in SRS
    
    # fit outcome model in SRS
    outcome = glm(y.s ~ T.s + propscore.srs, 
                  data = data.use.s,
                  family = stats::gaussian)
    coef.srs[i,] = outcome$coefficients
    
    # fit outcome model with true propensity score (weights only in outcome model)
    data.use1 = cbind(data.use, pi.b)
    outcome = svyglm(y.b ~ T.b + pi.b, 
                     design = survey::svydesign(ids = ~0, weights = data.use$htweight, data = data.use1),
                     family = stats::gaussian)
    coef.wt.outcome.trueprop[i,] = outcome$coefficients
    
    # fit outcome model with true propensity score (no weights)
    outcome = glm(y.b ~ T.b + pi.b, 
                  data = data.use1,
                  family = stats::gaussian)
    coef.nowt.trueprop[i,] = outcome$coefficients
    
    # fit outcome model in SRS with true propensity score
    data.use.s$pi.s = pi.s
    outcome = glm(y.s ~ T.s + pi.s, 
                  data = data.use.s,
                  family = stats::gaussian)
    coef.srs.trueprop[i,] = outcome$coefficients
    
  
  coef.check = rbind(apply(coef.srs,2,mean),
                     apply(coef.nowt, 2, mean),
                     apply(coef.wt.prop, 2, mean),
                     apply(coef.wt.outcome, 2, mean),
                     apply(coef.wt2, 2, mean),
                     apply(coef.nowt.trueprop, 2, mean),
                     apply(coef.wt.outcome.trueprop, 2, mean),
                     apply(coef.srs.trueprop, 2, mean))[,2]
  
  coef.check

  }
  
  coef.check = rbind(apply(coef.srs,2,mean),
                     apply(coef.nowt, 2, mean),
                     apply(coef.wt.prop, 2, mean),
                     apply(coef.wt.outcome, 2, mean),
                     apply(coef.wt2, 2, mean),
                     apply(coef.nowt.trueprop, 2, mean),
                     apply(coef.wt.outcome.trueprop, 2, mean),
                     apply(coef.srs.trueprop, 2, mean),
                     apply(coef.wt.outcome.Kprop, 2, mean),
                     apply(coef.wt2.Kprop, 2, mean),
                     apply(coef.wt.outcome.Kwt, 2, mean),
                     apply(coef.wt2.Kwt, 2, mean))[,2]
  
  coef.check
}

sfInit(parallel = TRUE,cpus = 40) # edit number of CPUS based on availability
sfExportAll()

sfLibrary(tidyverse)
sfLibrary(survey)
sfLibrary(splines)
sfLibrary(rlecuyer)
sfLibrary(rje)
sfClusterSetupRNG(type="RNGstream", seed = 1191913175)

start = Sys.time()

simresults = sfLapply(beta.list,sim)

end = Sys.time()
end-start

save(simresults, file = paste0(path,"coef",scenario,".RData"))


sfStop()


