#### How to estimate sampling weights using NHANES
#### Olivia Bernstein

library(tidyverse)
library(survey)
library(Hmisc)
library(estweight) # install with devtools::install_github("oliviabern/estweight")

## function from Dan Gillen's 211 code for estimating the robust sandwich se estimate
robust.se.glm <- function(glm.obj){
  ## 	Compute robust (sandwich) variance estimate
  robust.cov <- robust.vcov.glm(glm.obj)
  
  ##	Format the model output with p-values and CIs
  s <- summary( glm.obj) 
  robust.se <- sqrt( diag( robust.cov )) 
  z <- glm.obj$coefficients/robust.se
  p <- 2*pnorm( -abs( z ) ) 
  ci95.lo <- glm.obj$coefficients - qnorm( .975 ) * robust.se
  ci95.hi <- glm.obj$coefficients + qnorm( .975 ) * robust.se
  rslt <- cbind( glm.obj$coefficients, robust.se, ci95.lo, ci95.hi, z, p ) 
  dimnames(rslt)[[2]] <- c( dimnames( s$coefficients )[[2]][1], "Robust SE", "ci95.lo", "ci95.hi", "z value", "Pr(>|z|)" ) 
  rslt 
}
robust.vcov.glm <- function(glm.obj){
  if (is.matrix(glm.obj$x)) 
    xmat<-glm.obj$x
  else {
    mf<-model.frame(glm.obj)
    xmat<-model.matrix(terms(glm.obj),mf)		
  }
  umat <- residuals(glm.obj,"working")*glm.obj$weights*xmat
  modelv<-summary(glm.obj)$cov.unscaled
  robust.cov <- modelv%*%(t(umat)%*%umat)%*%modelv
  dimnames( robust.cov ) <- dimnames( vcov(glm.obj) )
  return( robust.cov )
}

#### Estimate sampling weights ####

# load nhanes data (code for creating nhanes_rep is in my estweight_simulation study repository)
nhanes_rep_full = read.csv("~/nhanes_rep.csv") # NHANES-REP

# load cleaned NACC data (from applied_example_clean_NACC.R)
load("~/nacc.RData")

# combine nursing (n=4) and assisted living (n=55) to allow for imputation
nacc$living = ifelse(nacc$living=="nursing","assistedliving",nacc$living)
# combine 4/completely dependence (n=7) and 3/requires some indpenedence (n=87)
nacc$independence = ifelse(nacc$independence==4,3,nacc$independence)

# prep nhanes_rep data to combine with nacc and create factors
# remove people under 60 from nhanes_rep
nhanes_rep = subset(nhanes_rep_full, age >= 60) %>% 
  dplyr::select(-kidney, -liver, -cancer, -exercise, -sleep) 

nhanes_rep$majordep = as.integer(nhanes_rep$majordep)#-1 # DOUBLE CHECK THIS
nh.factors = colnames(nhanes_rep)[c(3:10)]
nhanes_rep[nh.factors] <- lapply(nhanes_rep[nh.factors] , factor)

# compare demographics of NACC and NHANES (60+)
summary(nacc$age)
summary(nhanes_rep$age)

prop.table(table(nacc$race_eth)) # 74% NH White
prop.table(table(nhanes_rep$race_eth)) # 75% NH White

prop.table(table(nacc$edu, useNA = 'ifany')) # 63% college educated
prop.table(table(nhanes_rep$edu, useNA = 'ifany')) # 29% college educated

prop.table(table(nacc$female, useNA = 'ifany')) # 65% female
prop.table(table(nhanes_rep$female, useNA = 'ifany')) # 54% female

### Impute missing values in NACC data ####

set.seed(78951578)

# remove subjects who are missing the response faq.total
# removes 809 subjects for a total of 14,358 subjects
nacc.toimpute = nacc[!is.na(nacc$faq.total),]

factors = colnames(nacc.toimpute)[c(5:23)]
nacc.toimpute[factors] <- lapply(nacc.toimpute[factors] , factor)

# impute datasets
f = aregImpute(as.formula(paste0("~",paste0(colnames(nacc.toimpute)[c(2,4:23)],collapse = "+"))),
               data = nacc.toimpute, n.impute = 5)

M = f$n.impute
coef.save <- coef.save.onlyprop <- coef.save.onlyout <- coef.save.non <- matrix(NA, nrow = M, ncol = 3)
se.prop.save <- se.db.save <- se.naive.save <- matrix(NA, nrow = M, ncol = 3) # for 3 coefficients
se.db.save.onlyout <- se.sand.save.onlyprop <- se.sand.save.non <- matrix(NA, nrow = M, ncol = 3)

#### Run analysis on each imputed dataset #####

for(m in 1:M){
  
  print(m)
  
  #### Prep data ####
  
  # the m-th imputed dataset
  imputed = impute.transcan(f, imputation=m, data=nacc.toimpute, 
                                list.out=TRUE, pr=FALSE, check=FALSE)
  
  completed = nacc.toimpute
  completed[names(imputed)] <- imputed
  
  ## get dataset to use for estimating sampling weights 
  nacc.wts = completed %>% dplyr::select(SEQN, age_matching, female, edu,
                                           race_eth, highBP, diabetes, 
                                           CHF, dep, prescription) %>% 
    rename(age = age_matching, 
           majordep = dep)
  
  ## get dataset for analysis
  nacc.analysis1 = completed %>% dplyr::select(SEQN, age, female, edu,
                                                race_eth, highBP, diabetes, 
                                                CHF, dep, prescription,
                                                living, SPtype, married, independence,
                                                familyhx_cogimpair, tobac100,
                                                stroke, heartattack, seizures, thyroid,
                                                vitaminE, faq.total)
  #### Estimate sampling weights ####
  
  # create indicator for biased sample membership
  nacc.wts$biased = 1
  nhanes_rep$biased = 0
  
  # combine 2 datasets
  comb = rbind(nacc.wts, nhanes_rep)
  
  # convert biased to a factor
  comb$biased = as.factor(comb$biased)
  
  # remove spaces in factor names
  comb$race_eth = plyr::revalue(comb$race_eth, c("NH White" = "NHWhite", "NH Asian" = "NHAsian", "NH Black" = "NHBlack"))
  comb$edu = plyr::revalue(comb$edu, c("less than 12" = "lessthan12", "college grad or higher" = "collegeplus", 
                                       "high school/GED" = "highschool", "some college" = "somecollege"))
  
  # create smallest and largest formula statements
  preds = colnames(comb)[c(-1,-ncol(comb))]
  resp  = colnames(comb)[ncol(comb)]
  
  form.min = as.formula(paste(resp," ~ 1"))
  form.max = as.formula(paste(resp,"~",paste(preds,collapse = "+"),
                              "+ I(age^2)")) # add second order terms
  
  # form.max = as.formula(paste(resp,"~(",paste(preds,collapse = "+"),
  #                             ")^2",
  #                             "+ I(age^2)")) # if you want interaction terms
  
  # do stepwise AIC (forward selection)
  fit.min = glm(formula = form.min, family = binomial, data = comb)
  forward = step(fit.min,scope=list(lower=form.min,upper=form.max), 
                 direction="forward", trace = 0)
  
  
  # fit selected model
  estwt.form = formula(forward)
  estwt.fit = glm(formula = estwt.form, family = binomial, data = comb)
  
  # grab estimated HT weights
  prob.bias = fitted(estwt.fit, "response")       # est. probability biased == 1
  htweight_unnorm = (1-prob.bias)/prob.bias       # HT weight
  htweight = (htweight_unnorm/sum(htweight_unnorm))[comb$biased==1] # standardize to sum to 1
  
  estweights = data.frame(SEQN = comb$SEQN[comb$biased==1],
                          htweight = htweight)
  nacc.analysis = merge(nacc.analysis1,estweights, by = "SEQN")

  #### Estimate propensity score ####
  cont_vars_prop = colnames(nacc.analysis)[2] 
  fact_vars_prop = colnames(nacc.analysis)[3:20] 
  
  form_min_prop = stats::as.formula(paste(colnames(nacc.analysis)[21],"~ 1"))
  
  interactions = paste0(
    "(",paste(c(cont_vars_prop, fact_vars_prop),collapse = "+"),
    ")^2")
  form_max_prop = stats::as.formula(paste0(colnames(nacc.analysis)[21],"~",
                                           paste0(interactions, "+",
                                                  paste0(paste0("I(",cont_vars_prop,"^2)"),collapse = "+"))))
  
  fit_min = survey::svyglm(form_min_prop, 
                           design = survey::svydesign(ids = ~0, 
                                                      weights = nacc.analysis$htweight, 
                                                      data = nacc.analysis),
                           family = stats::quasibinomial)
  attach(nacc.analysis)
  forward = stats::step(fit_min,scope=list(lower=form_min_prop,upper=form_max_prop),
                        direction="forward", trace = 0)
  detach(nacc.analysis)
  
  # fit selected model
  estprop_form = stats::formula(forward)
  estprop = survey::svyglm(estprop_form, 
                           design = survey::svydesign(ids = ~0, weights = nacc.analysis$htweight, data = nacc.analysis),
                           family = stats::quasibinomial)
  
  nacc.analysis$propscore = fitted(estprop, "response")
  
  #### fit outcome model ####
  outcome_family = stats::gaussian
  outcome = svyglm(faq.total ~ vitaminE + propscore, 
                   design = survey::svydesign(ids = ~0, weights = nacc.analysis$htweight, data = nacc.analysis),
                   family = outcome_family)
  coef.save[m,] = outcome$coefficients
  
  # standard robust SE estimate
  se.db.save[m,] = summary(outcome)$coef[,2]

  # proposed SE estimate
  se.prop.save[m,] = treatmentSE(estwt_fit = estwt.fit, 
                               estprop_fit = estprop, 
                               fit_outcome = outcome, 
                               biased = comb$biased,
                               outcome_family = outcome_family)
  
  ## For comparison, could use the function from the estweight package
  # nhanes_rep_use = nhanes_rep %>%
  #   rename(dep = majordep, age_matching = age)
  # 
  # 
  # convPS(convSamp = completed, repSamp = nhanes_rep_use,
  #        sampwt_vars = c("age_matching", "female", "edu",
  #                       "race_eth", "highBP", "diabetes",
  #                       "CHF", "dep", "prescription"),
  #        PS_vars = c("age", "female", "edu",
  #                    "race_eth", "highBP", "diabetes",
  #                    "CHF", "dep", "prescription",
  #                    "living", "SPtype", "married", "independence",
  #                    "familyhx_cogimpair", "tobac100",
  #                    "stroke", "heartattack", "seizures", "thyroid"),
  #       treatment_var = "vitaminE", response_var = "faq.total")
  
  # naive (non sandwich) SE estimate
  outcome_naive = glm(faq.total ~ vitaminE + propscore, 
                      weights = nacc.analysis$htweight,
                      data = nacc.analysis,
                      family = stats::gaussian)
  se.naive.save[m,] = summary(outcome_naive)$coef[,2]
  
  #### fit outcome model (only prop weighted) ####
  outcome = glm(faq.total ~ vitaminE + propscore, 
                   data = nacc.analysis,
                   family = outcome_family)
  coef.save.onlyprop[m,] = outcome$coefficients
  se.sand.save.onlyprop[m,] = robust.se.glm(outcome)[,2]
  
  #### unweighted prop model ####
  fit_min = stats::glm(form_min_prop, 
                          data = nacc.analysis,
                           family = stats::binomial)
  
  
  attach(nacc.analysis)
  forward = stats::step(fit_min,scope=list(lower=form_min_prop,upper=form_max_prop),
                        direction="forward", trace = 0)
  detach(nacc.analysis)
  
  # fit selected model
  estprop_form = stats::formula(forward)
  estprop = stats::glm(estprop_form, data = nacc.analysis,
                           family = stats::binomial)
  
  nacc.analysis$propscore_unwt = fitted(estprop, "response")
  
  #### fit outcome model (only outcome weighted) ####
  outcome = svyglm(faq.total ~ vitaminE + propscore_unwt, 
                   design = survey::svydesign(ids = ~0, weights = nacc.analysis$htweight, data = nacc.analysis),
                   family = outcome_family)
  coef.save.onlyout[m,] = outcome$coefficients
  se.db.save.onlyout[m,] = summary(outcome)$coef[,2]
  
  #### fit outcome (no weighting) ####
  outcome = glm(faq.total ~ vitaminE + propscore_unwt, 
                data = nacc.analysis,
                family = outcome_family)
  coef.save.non[m,] = outcome$coefficients
  se.sand.save.non[m,] = robust.se.glm(outcome)[,2]
  
}


coef <- apply(coef.save, 2, mean)

var.db <- apply(se.db.save^2, 2, mean) + ((M+1)/M)*apply(coef.save, 2, var)
se.db <- sqrt(var.db)

var.prop <- apply(se.prop.save^2, 2, mean) + ((M+1)/M)*apply(coef.save, 2, var)
se.prop <- sqrt(var.prop)

var.naive <- apply(se.naive.save^2, 2, mean) + ((M+1)/M)*apply(coef.save, 2, var)
se.naive <- sqrt(var.naive)

ci = paste0(round(coef[2],2),
             " (",round(coef[2] - qnorm(.975)*se.prop[2],2)
             ,", ",round(coef[2] + qnorm(.975)*se.prop[2],2),")")

coef.onlyout <- apply(coef.save.onlyout, 2, mean)
se.onlyout <- sqrt(apply(se.db.save.onlyout^2, 2, mean) + ((M+1)/M)*apply(coef.save.onlyout, 2, var))
ci.onlyout = paste0(round(coef.onlyout[2],2),
                    " (",round(coef.onlyout[2] - qnorm(.975)*se.onlyout[2],2)
                    ,", ",round(coef.onlyout[2] + qnorm(.975)*se.onlyout[2],2),")")

coef.onlyprop <- apply(coef.save.onlyprop, 2, mean)
se.onlyprop <- sqrt(apply(se.sand.save.onlyprop^2, 2, mean) + ((M+1)/M)*apply(coef.save.onlyprop, 2, var))
ci.onlyprop = paste0(round(coef.onlyprop[2],2),
                    " (",round(coef.onlyprop[2] - qnorm(.975)*se.onlyprop[2],2)
                    ,", ",round(coef.onlyprop[2] + qnorm(.975)*se.onlyprop[2],2),")")

coef.non <- apply(coef.save.non, 2, mean)
se.non <- sqrt(apply(se.sand.save.non^2, 2, mean) + ((M+1)/M)*apply(coef.save.non, 2, var))
ci.non = paste0(round(coef.non[2],2),
                     " (",round(coef.non[2] - qnorm(.975)*se.non[2],2)
                     ,", ",round(coef.non[2] + qnorm(.975)*se.non[2],2),")")

results = matrix(cbind(ci.non, ci.onlyprop, ci.onlyout, ci),ncol = 1)

pt.est = matrix(cbind(coef.non[2], coef.onlyprop[2], coef.onlyout[2], coef[2]),ncol = 1)
sd =  matrix(cbind(se.non[2], se.onlyprop[2],
                   se.onlyout[2], se.prop[2]),ncol = 1)

ci.hi = est = matrix(cbind(coef.non[2] + qnorm(.975)*se.non[2], coef.onlyprop[2] + qnorm(.975)*se.onlyprop[2], 
                           coef.onlyout[2] + qnorm(.975)*se.onlyout[2], coef[2] + qnorm(.975)*se.prop[2]),ncol = 1)

rownames(results) = c("No weighting",
                      "Only propensity model weighted",
                      "Only outcome model weighted",
                      "Both models weighted")
colnames(results) = "Mean Est (95% CI)"

library("kableExtra")

kable(results, format = 'latex', booktabs = TRUE)

t = matrix(round(rbind(se.prop, se.db, se.naive),2)[,2], ncol = 1)
rownames(t) = c("proposed", "design based", "naive")
colnames(t) = "SE"
kable(t, format = 'latex', booktabs = TRUE)


# Forest plot of results

se.bars.vert = function(ci.lo, ci.hi, y, col, lwd = 3){
  arrows(ci.lo, y, ci.hi, y, length=0.05, angle=90, code=3, col = col,lwd = lwd)
}

# calculate CI
ci.lo = pt.est - qnorm(.975)*sd
ci.hi = pt.est + qnorm(.975)*sd

# create string for table
string = paste0(rnd(pt.est,2),
                " (",
                rnd(ci.lo,2),
                ", ", 
                rnd(ci.hi,2),")")

# Pick spots on y-axis for plotting
y = length(pt.est):1

png(filename="~/Dropbox/UCI/Research/Prediction Assessment/weighted_propensity_scores/applied_example/forest_plot.png", 
    width=3500, height=1500, res = 450)

# change plotting parameters to leave room for text
# on the left side 
par(oma = c(0,18,0,0) + 0.1, # adjust the second entry in the oma command to shift left side margin
    mar = c(5,0,4,1) + 0.1, ps = 10)


plot(pt.est,y, 
     xlim = range(c(ci.lo,ci.hi)),
     ylim = range(y) + c(-1,1), 
     yaxt = 'n', cex = 1.5,
     pch = 18, ylab = "", xlab = "Estimated Difference (95% CI) in FAQ Scores For Individuals \nWho Consume Vitamin E Supplements VS. Those Who Do Not",
     main = "Impact of Sampling Weights \non Estimated Difference in FAQ Scores")
abline(v = 0, lty = 3, lwd = 2, col = 'grey')
se.bars.vert(ci.lo, ci.hi, y = y, col = 'black')

# add text 
mtext(c("Estimated Difference in FAQ Scores (95% CI)"), 2, line = .5, 
      las = 2, at = max(y)+1, adj = 1, font = 2)
mtext(string, 2, line = .5, # adjust the line value to shift text left and right
      las = 2, at = y, adj = 1)

mtext(c("No weights",
        "Weight prop. model",
        "Weight outcome model",
        "Weight both models"), 
      2, line = 8, font = 1,
      las = 2, at = y, adj = 1)

dev.off()


#### Table 1 to show sample bias and impact of weighting ####
# use 1 imputation
# don't want to include missing data because you need a complete data set for estimating weights
# otherwise there will be a different sample size for each model in step()

source("~/Dropbox/UCI/Research/Prediction Assessment/Estimate C2C weights/applied example/Weighted Table 1 Functions.R")
source("~/Dropbox/UCI/Research/AD_trial_dyads/Completer analysis/table1functions.R")

## Estimate weights ##

imputed1 = impute.transcan(f, imputation=1, data=nacc.toimpute, 
                          list.out=TRUE, pr=FALSE, check=FALSE)

completed1 = nacc.toimpute
completed1[names(imputed1)] <- imputed1

## get dataset to use for estimating sampling weights 
nacc.wts.tab1 = completed1 %>% dplyr::select(SEQN, age_matching, female, edu,
                                       race_eth, highBP, diabetes, 
                                       CHF, dep, prescription) %>% 
  rename(age = age_matching, 
         majordep = dep)

nhanes_rep.tab1 = nhanes_rep

# create indicator for biased sample membership
nacc.wts.tab1$biased = 1
nhanes_rep.tab1$biased = 0

# combine 2 datasets
comb = rbind(nacc.wts.tab1, nhanes_rep.tab1)

# convert biased to a factor
comb$biased = as.factor(comb$biased)

# remove spaces in factor names
comb$race_eth = plyr::revalue(comb$race_eth, c("NH White" = "NHWhite", "NH Asian" = "NHAsian", "NH Black" = "NHBlack"))
comb$edu = plyr::revalue(comb$edu, c("less than 12" = "lessthan12", "college grad or higher" = "collegeplus", 
                                     "high school/GED" = "highschool", "some college" = "somecollege"))

# create smallest and largest formula statements
preds = colnames(comb)[c(-1,-ncol(comb))]
resp  = colnames(comb)[ncol(comb)]

form.min = as.formula(paste(resp," ~ 1"))
form.max = as.formula(paste(resp,"~",paste(preds,collapse = "+"),
                            "+ I(age^2)")) # add second order terms

# do stepwise AIC (forward selection)
fit.min = glm(formula = form.min, family = binomial, data = comb)
forward = step(fit.min,scope=list(lower=form.min,upper=form.max), 
               direction="forward", trace = 0)


# fit selected model
estwt.form = formula(forward)
estwt.fit = glm(formula = estwt.form, family = binomial, data = comb)

# grab estimated HT weights
prob.bias = fitted(estwt.fit, "response")       # est. probability biased == 1
htweight_unnorm = (1-prob.bias)/prob.bias       # HT weight
htweight = htweight_unnorm/sum(htweight_unnorm) # standardize to sum to 1

htweight.nacc = htweight[comb$biased == 1]

## nhanes summary
nh.sum = rbind(wt.sum(nhanes_rep.tab1$age),
               wt.prop(nhanes_rep.tab1$edu == 'less than 12'),
               wt.prop(nhanes_rep.tab1$edu == 'high school/GED'),
               wt.prop(nhanes_rep.tab1$edu == 'some college'),
               wt.prop(nhanes_rep.tab1$edu == 'college grad or higher'),
               wt.prop(nhanes_rep.tab1$race_eth == 'NH White'),
               wt.prop(nhanes_rep.tab1$race_eth == 'Hispanic'),
               wt.prop(nhanes_rep.tab1$race_eth == 'NH Asian'),
               wt.prop(nhanes_rep.tab1$race_eth == 'NH Black'),
               wt.prop(nhanes_rep.tab1$race_eth == "Other"),
               wt.prop(nhanes_rep.tab1$female == 1),
               wt.prop(nhanes_rep.tab1$highBP == 1),
               wt.prop(nhanes_rep.tab1$diabetes == 1),
               wt.prop(nhanes_rep.tab1$CHF == 1),
               wt.prop(nhanes_rep.tab1$majordep == 1),
               wt.prop(nhanes_rep.tab1$prescription == 1))


nacc.sum = rbind(wt.sum.meandiff(nacc.wts.tab1$age,nhanes_rep.tab1$age),
               wt.prop.meandiff(nacc.wts.tab1$edu == 'less than 12',nhanes_rep.tab1$edu == 'less than 12'),
               wt.prop.meandiff(nacc.wts.tab1$edu == 'high school/GED',nhanes_rep.tab1$edu == 'high school/GED'),
               wt.prop.meandiff(nacc.wts.tab1$edu == 'some college',nhanes_rep.tab1$edu == 'some college'),
               wt.prop.meandiff(nacc.wts.tab1$edu == 'college grad or higher',nhanes_rep.tab1$edu == 'college grad or higher'),
               wt.prop.meandiff(nacc.wts.tab1$race_eth == 'NH White',nhanes_rep.tab1$race_eth == 'NH White'),
               wt.prop.meandiff(nacc.wts.tab1$race_eth == 'Hispanic',nhanes_rep.tab1$race_eth == 'Hispanic'),
               wt.prop.meandiff(nacc.wts.tab1$race_eth == 'NH Asian',nhanes_rep.tab1$race_eth == 'NH Asian'),
               wt.prop.meandiff(nacc.wts.tab1$race_eth == 'NH Black',nhanes_rep.tab1$race_eth == 'NH Black'),
               wt.prop.meandiff(nacc.wts.tab1$race_eth == "Other",nhanes_rep.tab1$race_eth == "Other"),
               wt.prop.meandiff(nacc.wts.tab1$female == 1,nhanes_rep.tab1$female == 1),
               wt.prop.meandiff(nacc.wts.tab1$highBP == 1,nhanes_rep.tab1$highBP == 1),
               wt.prop.meandiff(nacc.wts.tab1$diabetes == 1,nhanes_rep.tab1$diabetes == 1),
               wt.prop.meandiff(nacc.wts.tab1$CHF == 1,nhanes_rep.tab1$CHF == 1),
               wt.prop.meandiff(nacc.wts.tab1$majordep == 1,nhanes_rep.tab1$majordep == 1),
               wt.prop.meandiff(nacc.wts.tab1$prescription == 1,nhanes_rep.tab1$prescription == 1))

nacc.sum.wt = rbind(wt.sum.meandiff(nacc.wts.tab1$age,nhanes_rep.tab1$age, htweight.nacc),
                    wt.prop.meandiff(nacc.wts.tab1$edu == 'less than 12',nhanes_rep.tab1$edu == 'less than 12', htweight.nacc),
                    wt.prop.meandiff(nacc.wts.tab1$edu == 'high school/GED',nhanes_rep.tab1$edu == 'high school/GED', htweight.nacc),
                    wt.prop.meandiff(nacc.wts.tab1$edu == 'some college',nhanes_rep.tab1$edu == 'some college', htweight.nacc),
                    wt.prop.meandiff(nacc.wts.tab1$edu == 'college grad or higher',nhanes_rep.tab1$edu == 'college grad or higher', htweight.nacc),
                    wt.prop.meandiff(nacc.wts.tab1$race_eth == 'NH White',nhanes_rep.tab1$race_eth == 'NH White', htweight.nacc),
                    wt.prop.meandiff(nacc.wts.tab1$race_eth == 'Hispanic',nhanes_rep.tab1$race_eth == 'Hispanic', htweight.nacc),
                    wt.prop.meandiff(nacc.wts.tab1$race_eth == 'NH Asian',nhanes_rep.tab1$race_eth == 'NH Asian', htweight.nacc),
                    wt.prop.meandiff(nacc.wts.tab1$race_eth == 'NH Black',nhanes_rep.tab1$race_eth == 'NH Black', htweight.nacc),
                    wt.prop.meandiff(nacc.wts.tab1$race_eth == "Other",nhanes_rep.tab1$race_eth == "Other", htweight.nacc),
                    wt.prop.meandiff(nacc.wts.tab1$female == 1,nhanes_rep.tab1$female == 1, htweight.nacc),
                    wt.prop.meandiff(nacc.wts.tab1$highBP == 1,nhanes_rep.tab1$highBP == 1, htweight.nacc),
                    wt.prop.meandiff(nacc.wts.tab1$diabetes == 1,nhanes_rep.tab1$diabetes == 1, htweight.nacc),
                    wt.prop.meandiff(nacc.wts.tab1$CHF == 1,nhanes_rep.tab1$CHF == 1, htweight.nacc),
                    wt.prop.meandiff(nacc.wts.tab1$majordep == 1,nhanes_rep.tab1$majordep == 1, htweight.nacc),
                    wt.prop.meandiff(nacc.wts.tab1$prescription == 1,nhanes_rep.tab1$prescription == 1, htweight.nacc))

out.sum = cbind(nh.sum, nacc.sum, nacc.sum.wt)
rownames(out.sum) = c("Age",
                      paste0("Educ: ", c("Less than high school", "High school", "Some college", "College or higher")),
                      paste0("Race/Ethnicity: ", c("NH White", "Hispanic", "NH Asian", "NH Black", "Other")),"Female",
                      "High blood pressure", "Diabetes", "CHF", 
                      "Major depression", "Prescription meds")
colnames(out.sum) = c("Unweighted","Unweighted","Weighted")

# kable.booktab(out.sum,nvaringroup = c(1,4,5,6),
#               fittopagewidth = TRUE, escape = TRUE)

k = kable(out.sum, format = 'latex', booktabs = TRUE,
          linesep = get.linesep(c(1,4,5,6)),
          caption = "ADD")%>%
  kable_styling(latex_options = c("hold_position", "scale_down"))

add_header_above(k, c(" "=1,"NHANES-REP" = 1,
                      "NACC" = 2))
