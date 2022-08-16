#### How to estimate sampling weights using NHANES
#### Olivia Bernstein Morgan

library(tidyverse)

#### Estimate sampling weights ####

# load NACC data
nacc_full = read.csv("~/investigator_nacc55.csv")

# number of ADRCs in dataset
length(unique(nacc_full$NACCADC)) # 41 different ADRCs

nacc = nacc_full %>% dplyr::select(NACCID, NACCVNUM,
                                   #VISITMO, VISITDAY, 
                                   VISITYR,
                                   NACCAGE, 
                                   SEX, HISPANIC, RACE,
                                   INRELTO,
                                   EDUC, MARISTAT, NACCLIVS,
                                   INDEPEND, RESIDENC,
                                   NACCFAM,
                                   #NACCFADM, NACCAM,
                                   NACCUDSD,
                                   ANYMEDS, 
                                   starts_with("DRUG"), # here
                                   TOBAC100, 
                                   CVCHF, CONGHRT,
                                   CBSTROKE,
                                   CVHATT, MYOINF,
                                   SEIZURES, 
                                   DIABETES, DIABET,
                                   HYPERTEN, HXHYPER,
                                   DEP,
                                   THYDIS, THYROID, 
                                   BILLS, TAXES, SHOPPING, GAMES,
                                   STOVE, MEALPREP, EVENTS,
                                   PAYATTN, REMDATES, TRAVEL) %>% 
  filter(NACCVNUM == 1,
         NACCAGE >= 60, 
         NACCUDSD == 1) %>%
  mutate(age_matching = ifelse(NACCAGE > 80,80,NACCAGE),
         female = ifelse(SEX==2,1,0),
         race_eth = ifelse(HISPANIC==1,"Hispanic", 
                           ifelse(RACE==1, "NH White",
                                  ifelse(RACE==2,"NH Black",
                                         ifelse(RACE==5,"NH Asian",
                                                ifelse(RACE==99,NA,"Other"))))),
         edu = ifelse(EDUC==99, NA,
                      ifelse(EDUC < 12, "less than 12",
                             ifelse(EDUC==12, "high school/GED",
                                    ifelse(EDUC < 16, "some college","college grad or higher")))),
         living = ifelse(RESIDENC == 1,"privateres",
                         ifelse(RESIDENC == 2, "retirement",
                                ifelse(RESIDENC == 3, "assistedliving",
                                       ifelse(RESIDENC == 4, "nursing","other_missing")))),
         diabetes = ifelse((DIABETES %in% 1:2)|(DIABET %in% 1:3),1,
                           ifelse((DIABETES == 0) | (DIABET == 0),0,NA)),
         CHF = ifelse((CVCHF %in% 1:2)|CONGHRT==1,1,
                      ifelse(CVCHF == 0 | CONGHRT == 0,0,NA)),
         highBP = ifelse((HYPERTEN %in% 1:2) | HXHYPER == 1, 1,
                         ifelse(HYPERTEN == 0 | HXHYPER == 0, 0, NA)),
         SPtype = ifelse(INRELTO==-4,NA,INRELTO),
         independence = ifelse(INDEPEND == 9, NA, INDEPEND),
         familyhx_cogimpair = ifelse(NACCFAM == -4, NA, NACCFAM), 
         tobac100 = ifelse(TOBAC100==9,NA,TOBAC100), 
         stroke = ifelse(CBSTROKE==9,NA,CBSTROKE),
         heartattack = ifelse((CVHATT %in% 1:2)|(MYOINF == 1), 1, 
                              ifelse((CVHATT == 0)|(MYOINF == 0),0,NA)),
         seizures = ifelse(SEIZURES == 9, NA, SEIZURES), 
         thyroid = ifelse((THYDIS == 1)|(THYROID %in% 1:2),1,
                          ifelse((THYDIS == 0)|(THYROID == 0),0,NA))) %>% 
  rename(SEQN = NACCID, 
         age = NACCAGE,
         dep = DEP, 
         married = MARISTAT) %>% 
  dplyr::select(SEQN, VISITYR,
                age_matching, age, female, race_eth,
                edu, living,
                diabetes, highBP,
                dep, CHF,
                ANYMEDS, starts_with("DRUG"),
                BILLS, TAXES, SHOPPING, GAMES,
                STOVE, MEALPREP, EVENTS,
                PAYATTN, REMDATES, TRAVEL,
                SPtype, married, independence,
                familyhx_cogimpair, tobac100, stroke,
                  heartattack, seizures, 
                thyroid)

# # check family history of dominant mutation & mutation
# dominant_mutation = as.factor(nacc$NACCFADM)
# levels(dominant_mutation) = c("no/unknown","yes")
# any_mutation = as.factor(nacc$NACCAM)
# levels(any_mutation) = c("Not available", "No", "Yes, APP", "Yes, PS-1", "Yes, PS-2", "Yes, other","Unknown")
# 
# table(dominant_mutation, any_mutation)
# table(any_mutation, dominant_mutation)


# NACC includes OTC drugs, we want an indicator of which patients take prescription drugs
# we will remove common OTC drugs from the drug list
# we will start with the top 10 OTC drugs reported in NHANES 

a1 = nacc %>% dplyr::select(starts_with("DRUG"))
a2 = apply(a1, 2, as.character)
alldrugs = array(a2, dim = c(nrow = nrow(a2)*ncol(a2),1))
uniqdrugs = unique(alldrugs)
# 1352 unique drugs in NACC (for people 60+ who are cognitively normal, and )

t = data.frame(table(alldrugs))
t[order(t[,"Freq"],decreasing = TRUE),][1:40,]
# checked the top 40 drugs in NACC
# it includes several vitamins (including ascorbic acid so I removed that as well)

# 10 most popular OTC drugs in NHANES
otc_NHANES =  c("IBUPROFEN", "FLUTICASONE", "CETIRIZINE", "POTASSIUM CHLORIDE", 
         "ESOMEPRAZOLE", "NAPROXEN", "ASPIRIN", "TRIAMCINOLONE","ADVIL",
         "ZYRTEC")
# other common pain killers + vitamins 
otc = c(otc_NHANES, "TYLENOL", "ACETAMINOPHEN", "IBUPROFEN",
              "MIDOL","MOTRIN", "ACETAMINOPHEN", "LORATADINE","CLARITIN",
              "VITAMIN", "ASCORBIC ACID")


# check which NACC drugs contain the terms in "otc"
torem = list()
for(i in 1:length(otc)){
  torem[[i]] = as.matrix(uniqdrugs[which(grepl(otc[i],uniqdrugs,ignore.case = TRUE))], ncol = 1)
}

torem = do.call("rbind",torem)

# there are some combination drugs that include both OTC and prescription ingredients
# for example, "ACETAMINOPHEN-CODEINE"
# these drugs should still be considered prescription and stay in the drugs list
# I looked through all of the OTC drugs in "torem" and checked for prescription drugs
# the list of prescription drugs I found is below in "presc_w_otc"

presc_w_otc = c("HYDROCODONE", "SUMATRIPTAN", "BUTALBITAL",
            "CARISOPRODOL", "METHOCARBAMOL", "DIPYRIDAMOLE", "CODEINE",
            "ORPHENADRINE", "NYSTATIN", "OXYCODONE", "TRAMADOL", "PROPOXYPHENE",
            "DIHYDROCODEINE", "DESLORATADINE")

# otc/presc combos to take out of list
combos = matrix(NA, nrow = nrow(torem), ncol = length(presc_w_otc))
for(i in 1:length(presc_w_otc)){
  combos[,i] = grepl(presc_w_otc[i], torem, ignore.case = TRUE)
}
combos.total = apply(combos,1,max)==1

# the final list of OTC drugs to remove (without otc/prescription combos)
torem.final = unique(torem[!combos.total]) # 59 drugs to remove

# list of the remaining prescription drugs
presc1 = uniqdrugs[!(uniqdrugs %in% torem.final)]
presc = presc1[-1] # remove blank entry "" from list

## create prescription indicator
# count of prescription drugs
np1 = nacc %>% dplyr::select(SEQN, starts_with("DRUG"))
np = apply(np1, 1,function(x){sum(x %in% presc)})
anypresc = as.integer(np>0)

# treat people with ANYMEDS == -4 as NA
prescription = ifelse(nacc$ANYMEDS == -4, NA, anypresc)

# add into nacc data frame
nacc$prescription = prescription

## get vitamin E
uniqdrugs[which(grepl("VITAMIN E",uniqdrugs,ignore.case = TRUE))]
nE = apply(np1, 1,function(x){sum(x == "VITAMIN E")})
nacc$vitaminE = ifelse(nacc$ANYMEDS == -4, NA, as.integer(nE>0))


# GET FAQ (functional activities questionnaire)
# if answer is "not applicable", I will give them a score of 0
# 
faq = nacc %>% dplyr::select(BILLS, TAXES, SHOPPING, GAMES,
                             STOVE, MEALPREP, EVENTS,
                             PAYATTN, REMDATES, TRAVEL) %>% 
  mutate(across(everything(),~ifelse(.x %in% c(-4, 9), NA, .x)),
         across(everything(),~ifelse(.x == 8, NA, .x)) # can replace with 0
         )
# apply(faq,2,function(x){sum(x==8, na.rm = TRUE)})

# get total scores with imputation
faq.nNA = apply(faq, 1, function(x){sum(is.na(x))})
faq.total.na.rm = apply(faq, 1, function(x){sum(x, na.rm = TRUE)})

faq.total = ifelse(faq.nNA == 0, faq.total.na.rm,
                   ifelse(faq.nNA <= 3, faq.total.na.rm*10/(10-faq.nNA),NA))

# faq[which((faq.total==5)&(faq.nNA==2)),]
# table(faq.total, faq.nNA, useNA = 'ifany')

nacc$faq.total = faq.total

# remove unnecessary covariates
nacc = nacc %>% dplyr::select(-ANYMEDS, -starts_with("DRUG"),
                                         -BILLS, -TAXES, -SHOPPING,
                                         -GAMES, -STOVE, -MEALPREP, -EVENTS,
                                         -PAYATTN, -REMDATES, -TRAVEL)


save(nacc, file = "~/nacc.RData")

