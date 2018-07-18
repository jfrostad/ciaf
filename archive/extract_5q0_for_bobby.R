## ///////////////////////////////////////////////// ##
## ROY BURSTEIN
## DECEMBER 2017
## extract 5q0 by year-survey-admin1 for bobby
## ///////////////////////////////////////////////// ##

## ///////////////////////////////////////////////// ##
## SETUP

## load packages
library(data.table)
library(survival)
library(Biograph)
library(rms)

## set your root directory
root <- 'C:/Users/royburst/Desktop/cbh_cgf'

# load data
d <- readRDS(sprintf('%s/CBH_CGF_extract.RDS',root))
d <- setDT(d)


## ///////////////////////////////////////////////// ##
## Age bin cutpoints to use for synthetic cohort tabulation
ab_times <- data.frame(
  tstart = c(0,29/30,6,12,24,36,48),
  tstop  = c(29/30,6,12,24,36,48,60),
  ab     = c('nn','pnn1','pnn2','1y','2y','3y','4y') )
ab_times$ab <- paste0(1:length(ab_times$ab),"_",ab_times$ab)
ab_times$tmid <- ab_times$tstart+(ab_times$tstop-ab_times$tstart)/2



## ///////////////////////////////////////////////// ##
## Prep the data to be reshaped by age bin
d[,childage := interview_date_cmc-child_dob_cmc]  
d$childage[d$child_alive==0] <- d$child_age_at_death_months[d$child_alive==0]
d[,childage := childage + 0.0001]
d[,died  := child_alive==0]
d[,yrborn := child_dob_cmc/12+1900]

# only keep children born in the past 10 years 
# (so we can get 5q0 for last 5). Mostly reduce the data to quicken up the reshape
d <- subset(d, childage <= 120)



## ///////////////////////////////////////////////// ##
## reshape data by agebin so we
d <- d[,c('nid','source','year','yrborn','weight','cluster_number','strata','admin_1','childage','died'),with=FALSE]
dl <- data.table(survSplit(Surv(childage, died) ~.,cut=ab_times$tstart,data=d))
dl <- merge(dl,ab_times,by='tstart') # add on age bin information
dl[, year_ent := floor(yrborn + tstart/12)] # identify the year entering each bin
dl[,N:=1]

# keep only children-bin rows if entered age bin in the past five years
dl <- subset(dl, year-year_ent <=5) 



## ///////////////////////////////////////////////// ##
## collapse admin 1 by year and age
# get weights first
dl <- dl[,pw := weight/sum(weight), by = .(nid,ab,admin_1)]
# agebin specific probability of death
cbh <- dl[, .(q        = sum((died/N)*pw),         # probability of death before age 5
              died     = sum(died),                # total died
              childyrs = sum((tstop-tstart)/12)),  # Sample size
          by = .(nid,source,ab,admin_1)]           # aggregate to agebin-admin1

# collapse age bins to get q5
cbh <- cbh[, .(q = 1-prod(1-q), died = sum(died), childyrs = sum(childyrs)), 
           by = .(nid,source,admin_1)]
cbh <- cbh[order(-nid,admin_1)]


## ///////////////////////////////////////////////// ##
## Save and have a good day
write.csv(cbh,sprintf('%s/u5m_for_bobby.csv',root))


