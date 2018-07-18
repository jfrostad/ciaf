# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: JF
# Date: 06/12/2018
# Purpose: Prepare CGF data for the CIAF project
# source("/home/j/temp/jfrostad/ciaf/code/prep.R", echo=T)
#***********************************************************************************************************************

# ----CONFIG------------------------------------------------------------------------------------------------------------
# clear memory
rm(list=ls())

# runtime configuration
if (Sys.info()["sysname"] == "Linux") {
  j_root <- "/home/j/"
  h_root <- "/homes/jfrostad/"
  arg <- commandArgs()[-(1:3)] # First args are for unix use only

  if (length(arg)==0) {
    # arg <- c("IND", #current project iteration
    #          "8", #output version
    #          1) #number of cores provided to multicore functions
  }

  package_lib    <- sprintf('%s_code/_lib/pkg',h_root)
  ## Load libraries and  MBG project functions.
  .libPaths(package_lib)

  # necessary to set this option in order to read in a non-english character shapefile on a linux system (cluster)
  Sys.setlocale(category = "LC_ALL", locale = "C")

} else {
  j_root <- "J:"
  h_root <- "H:"
  # arg <- c("IND", #current project iteration
  #          "4", #output version
  #          1) #number of cores provided to multicore functions
}

# load packages
pacman::p_load(Biograph, data.table, faraway, grid, ggplot2, magrittr, mgcv, nipnTK, parallel,
               reshape2, rms, sf, SDMTools, stringr, survival, tidyverse, viridis)

# set working directories
home.dir <- file.path(j_root, "/temp/jfrostad/ciaf/")
  setwd(home.dir)

#run in parallel?
cores <- 6

#toggles
prep.data <- TRUE
summary.plot <- FALSE
run.model <- TRUE
run.scrap <- F

# z score standard cutoff to determine categories
z.cutoff <- -3

#data quality drops toggles
who.flag.drops <- T
who.flag.pct.drops <- T
age.ratio.drops <- T
age.dps.drops <- T
mean.dps.drops <- F
sd.drops <- T

#variables
categories <- c('none', 'wasting', 'underweight', 'stunting',  'wasting_underweight', 'stunting_underweight', 'all')
indicators <- paste0('i_', categories)
prevs <- paste0('p_', categories)

temp <- expand.grid(x=c('haz', 'waz', 'whz'), y=c('haz', 'waz', 'whz')) %>% as.data.table
temp[, names := paste0(y, '_', x)]
corr.vars <- unique(temp$names)
#***********************************************************************************************************************

# ----IN/OUT------------------------------------------------------------------------------------------------------------
###Input###
#raw data
data.dir <- file.path(home.dir, 'data')
raw.dt <- file.path(data.dir, 'raw', 'raw_data_bobby.csv') %>% fread

#master shapefiles
shapefile.dir <- file.path(j_root, "DATA/SHAPE_FILES/GBD_geographies/master/GBD_2016/master/shapefiles")
shapefile.version <- "GBD2016_mapping_final"
borders <- st_read(shapefile.dir, layer = shapefile.version)

###Output###
graph.dir <- file.path(home.dir, 'graphs')
#***********************************************************************************************************************

# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
##function lib##
#general functions#
central.function.dir <- file.path(h_root, "_code/_lib/functions")
# this pulls the general misc helper functions
file.path(central.function.dir, "misc.R") %>% source
#shared functions#
shared.function.dir <- file.path(j_root,  "temp/central_comp/libraries/current/r")
file.path(shared.function.dir, 'get_location_metadata.R') %>% source
file.path(shared.function.dir, 'get_ids.R') %>% source
file.path(shared.function.dir, 'get_covariate_estimates.R') %>% source

#custom fx
invLogit <- function(x) 1/(1+exp(-x))

##mapping##
mapPoints <- function(map.dt, year_start, year_end, var, scale.dir=1) {

  plot <- ggplot(borders) +
    geom_sf() +
    geom_point(data=map.dt[year > year_start & year < year_end], aes_string(x='long', y='lat', color=var))+
    scale_color_distiller(var, palette = "YlOrBr", direction=scale.dir) +
    ggtitle(var, subtitle=paste0("Year = ", year_start, "-", year_end)) +
    theme_minimal()

  return(plot)

}

scatterVars <- function(plot.dt, year.start, year.end, x.var, y.var, color.var,
                        colors=region.palette) {

  message('working on x:', x.var, '/y:', y.var)

  plot <- ggplot(data=plot.dt[year > year.start & year < year.end],
                 aes_string(x=x.var, y=y.var, color=color.var)) +
    geom_point() +
    scale_colour_manual(values=colors) +
    ggtitle(paste0(x.var, " vs. ", y.var),
            subtitle=paste0("Year = ", year.start, "-", year.end)) +
    theme_minimal()

  return(plot)

}

annotation_custom2 <-
  function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data)
  {
    layer(data = data, stat = StatIdentity, position = PositionIdentity,
          geom = ggplot2:::GeomCustomAnn,
          inherit.aes = TRUE, params = list(grob = grob,
                                            xmin = xmin, xmax = xmax,
                                            ymin = ymin, ymax = ymax))
  }

gamPreds <- function(input.dt, model.fam,
                     year.start, year.end,
                     pred.var, response.var,
                     k.value, by.var=NULL, weight.var=NULL,
                     transformeR=plogis) {

  mod.dt <- input.dt[year>=year.start & year<=year.end] %>% copy #subset years
  setnames(mod.dt, pred.var, 'x') #rename the dependent variable to make prediction plotting easier
  weight.vector <- ifelse(rep(is.null(weight.var), nrow(mod.dt)),
                          NULL,
                          mod.dt[, get(weight.var)]) #setup weight vector if weight var has been passed

  message('working on predictor:', pred.var, '/response:', response.var, "...(", year.start, "/", year.end, ")")

  form <- paste0(response.var, " ~ s(x, by=as.factor(", by.var, "), k=", k.value, ")") %>% as.formula
  message("formula: ", form)

  #run model
  mod <- gam(form,
             family=model.fam,
             weights=weight.vector,
             data=mod.dt)

  #setup dt for predictions and then duplicate it across your by variable
  pred.dt <- data.table(x=seq(0,1, by=0.025))

  duplicatoR <- function(value, var, dt) {

    out <- copy(dt)
    out[, (var) := value]
    return(out)

  }

  pred.dt <- lapply(unique(mod.dt[, get(by.var)]), duplicatoR,
                    by.var, pred.dt) %>% rbindlist

  new.vars <- c('pred', 'se', 'lower', 'upper')

  #generate predictions and CI
  pred.dt[, (new.vars[1:2]) := predict(mod, pred.dt, family=model.fam, type='link', se.fit=T)]
  pred.dt[, (new.vars[3]) := pred - 1.96*se]
  pred.dt[, (new.vars[4]) := pred + 1.96*se]
  pred.dt[, (new.vars) := lapply(.SD, transformeR), .SDcols=new.vars]


  #save the prediction var for labelling
  pred.dt[, predictor := pred.var]

  return(pred.dt)

}

scatterPreds <- function(plot.dt, year.start, year.end, color.var,
                         colors=region.palette, this.title, this.subtitle) {

  plot <- ggplot(data=plot.dt, aes_string(x='x', y='pred', color=color.var)) +
    geom_line() +
    geom_ribbon(aes_string(ymin='lower', ymax='upper', fill=color.var), alpha=.1) +
    facet_wrap(~predictor)+
    scale_colour_manual(values=colors) +
    scale_fill_manual(values=colors) +
    ggtitle(paste0(this.title, " Year = ", year.start, "-", year.end),
            subtitle=this.subtitle) +
    theme_minimal()

  return(plot)

}

#***********************************************************************************************************************

# ---PREP DATA----------------------------------------------------------------------------------------------------------
##prep data and create summary statistics of CGF indicators##
#bring in location information
locs <- get_location_metadata(location_set_id = 9)

#bring in covs
sdi <- get_covariate_estimates(covariate_id = 881)
  sdi <- sdi[, .(location_id, year_id, mean_value)]
  setnames(sdi, c('year_id', 'mean_value'), c('year', 'sdi'))
haqi <- get_covariate_estimates(covariate_id = 1099)
  haqi <- haqi[, .(location_id, year_id, mean_value)]
  setnames(haqi, c('year_id', 'mean_value'), c('year', 'haqi'))

#prep the DT
setnames(raw.dt, 'V1', 'id')

#create a cluster ID (cluster_num doesnt seem to be unique except within NID)
raw.dt[, cluster_id := paste0(nid, '_', cluster_number)]

#rename some vars
setnames(raw.dt, c('longnum', 'latnum', 'country'), c('long', 'lat', 'ihme_loc_id'))

#merge on location labels
raw.dt <- merge(raw.dt, locs[, .(ihme_loc_id, location_id, location_name, region_name, super_region_name)], by='ihme_loc_id')

#merge on covs
raw.dt <- merge(raw.dt, sdi, by=c('location_id', 'year'))
raw.dt <- merge(raw.dt, haqi, by=c('location_id', 'year'))
setkey(raw.dt, id) #rekey on row ID

#also create a color palette for plotting by regions
region.palette <- RColorBrewer::brewer.pal(unique(raw.dt$region_name), 'Paired')
region.palette[11] <- "#B3B3B3" #replace yellow with gray
region.palette <- c(region.palette, "#E5C494") #add on a 13th color (brown)
#***********************************************************************************************************************

# ---CALC U5M-----------------------------------------------------------------------------------------------------------
##prep data to reshape and calculate U5M
#Age bin cutpoints to use for synthetic cohort tabulation
ab.times <- data.table(
  tstart = c(0,29/30,6,12,24,36,48),
  tstop  = c(29/30,6,12,24,36,48,60),
  ab     = c('nn','pnn1','pnn2','1y','2y','3y','4y'))

ab.times[, ab := paste0(1:length(ab), "_", ab)]
ab.times[, tmid := tstart + (tstop - tstart)/2]

## Prep the data to be reshaped by age bin
u5m.dt <- copy(raw.dt)
u5m.dt[, childage := interview_date_cmc-child_dob_cmc]
u5m.dt[child_alive==0, childage := child_age_at_death_months] #if child is dead, record age at death
u5m.dt[, childage := childage + 0.0001] #offset zeros?
u5m.dt[, died  := child_alive==0]
u5m.dt[, yrborn := child_dob_cmc/12+1900]

# only keep children born in the past 10 years
# (so we can get 5q0 for last 5). Mostly reduce the data to quicken up the reshape
u5m.dt <- u5m.dt[childage <= 120]

## reshape data by agebin so we
u5m.split <- u5m.dt[,c('nid','source','year','yrborn','weight','cluster_number','strata','admin_1','childage','died')]
u5m.split <- survSplit(Surv(childage, died) ~ ., cut=ab.times$tstart, data=u5m.split) %>% as.data.table
u5m.split <- merge(u5m.split, ab.times, by='tstart') # add on age bin information
u5m.split[, year_ent := floor(yrborn + tstart/12)] # identify the year entering each bin
u5m.split[, N := 1]

# keep only children-bin rows if entered age bin in the past five years
u5m.split <- u5m.split[(year-year_ent)<=5]

## collapse admin 1 by year and age
# get weights first
u5m.split[, pw := weight/sum(weight), by = .(nid, ab, admin_1)]
# agebin specific probability of death
u5m.split <- u5m.split[, .(q        = sum((died/N)*pw),         # probability of death before age 5
                           died     = sum(died),                # total died
                           #survived = sum(died==0),             # total survivors TODO, need to figure out a better way to do this as it doesnt work currently
                           childyrs = sum((tstop-tstart)/12)),  # Sample size
                       by = .(nid, source, ab, admin_1)]           # aggregate to agebin-admin1

# collapse age bins to get q5
u5m.split <- u5m.split[, .(u5m = 1-prod(1-q),
                           died = sum(died),
                           #survived = sum(survived),
                           childyrs = sum(childyrs)),
                       by = .(nid, source, admin_1)]
u5m.dt <- u5m.split[order(-nid, admin_1)]

#bin the u5m values for small multiples plots
u5m.dt[, u5m_bin := cut(u5m, quantile(u5m, seq(0, 1, by = 0.25)))]
#***********************************************************************************************************************

# ---CALC CGF-----------------------------------------------------------------------------------------------------------
#drop data that is missing from the anthro module
#this data is merged between birth registry and anthro, so many births are missing the anthro data
cgf.dt <- na.omit(raw.dt, cols=c('HAZ', 'WAZ', 'WHZ'))
og.group.count <- cgf.dt %>% uniqueN(by=c('nid', 'source', 'admin_1')) #save the original count of groups to compare after data cleaning

#number of children surveyed
cgf.dt[, pop := sum(N), by=.(nid, source, admin_1)]

#some data also needs to be excluded/dropped based on CGF criteria (#TODO, find out more about these criteria)
cgf.dt[, exclude := 0]
cgf.dt[(HAZ_exclude == 1 | WAZ_exclude == 1 | WHZ_exclude == 1), exclude := 1]
cgf.dt[, drop := 0]
cgf.dt[(HAZ_drop == 1 | WAZ_drop == 1 | WHZ_drop == 1), drop := 1]

#do some summarizing and exploration
if (summary.plot==T) {
  #do some mapping
  mapPoints(cgf.dt, 2010, 2016, 'exclude', -1)
  mapPoints(cgf.dt, 2010, 2016, 'drop', -1)
  #xplore data distributions
  ggplot(cgf.dt, aes(x=child_height, fill=region_name)) + geom_density(alpha=.5) + scale_fill_manual(values=region.palette) + theme_minimal()
  ggplot(cgf.dt, aes(x=child_weight, fill=region_name)) + geom_density(alpha=.5) + scale_fill_manual(values=region.palette) + theme_minimal()
  ggplot(cgf.dt, aes(x=age_month, fill=region_name)) + geom_density(alpha=.5) + scale_fill_manual(values=region.palette) + theme_minimal()

  #simulate the amount of admin1s that would be dropped based on various criteria for percent lost
  criteria.list <- seq(0, 1, by=.05)

  simCriteria <- function(x, dt, by.list) {

    num <- dt[percent_lost > x] %>% uniqueN(by=by.list)
    out <- data.table(criteria=x, percent_lost=num/uniqueN(dt, by=by.list))

  }

  drop.dt <- lapply(criteria.list, simCriteria, dt=cgf.dt, by.list=c('nid', 'source', 'admin_1')) %>% rbindlist
}

#use a combined count of points that need to be dropped or excluded as an indicator of data quality within the admin1
cgf.dt[, drop_count := (drop + exclude) %>% sum, by=.(nid, source, admin_1)]
cgf.dt[, percent_lost := drop_count/pop]

#drop any admin1 groups that have more than 10% of rows being dropped due to data quality issues
drop.criteria <- .1
message("there are ",
        cgf.dt[percent_lost > drop.criteria] %>% uniqueN(by=c('nid', 'source', 'admin_1')),
        " groups that have more than 10% data quality issues. \n see admin 1s below:")
cgf.dt[percent_lost > drop.criteria, paste0(ihme_loc_id, "->", admin_1) %>% table]

#remove the points
if(who.flag.drops==TRUE) cgf.dt <- cgf.dt[!cgf.dt[(exclude == 1 | drop == 1)]]
if(who.flag.pct.drops==TRUE) cgf.dt <- cgf.dt[!(percent_lost > drop.criteria)]

#cleanup
cgf.dt[, names(cgf.dt)[names(cgf.dt) %like% 'drop'] := NULL]
cgf.dt[, names(cgf.dt)[names(cgf.dt) %like% 'exclude'] := NULL]
cgf.dt[, names(cgf.dt)[names(cgf.dt) %like% '_b'] := NULL]

#also remove the albanian survey
#no other datapoints in europe and these children seem to be significantly older, taller, and heavier than all others
cgf.dt <- cgf.dt[ihme_loc_id != "ALB"]

#further exploration of data quality
#checks based on this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5496063/

#age ratio: ratio of children aged 6–29 months over the number of children aged 30–59 months
#this ratio should be around 1 if the sample is of good quality
cgf.dt[, age_ratio := sum(age_month>=6 & age_month<30)/sum(age_month>=30 & age_month<60), by=.(nid, source, admin_1)]

#drop any admin1 groups that are >5 ratio
drop.criteria <- 4
message("there are ",
        cgf.dt[abs(age_ratio-1)>drop.criteria] %>% uniqueN(by=c('nid', 'source', 'admin_1')),
        " groups that have unreasonable age ratio. \n see admin 1s below:")
cgf.dt[abs(age_ratio-1)>drop.criteria, paste0(ihme_loc_id, "->", admin_1) %>% table]
if(age.ratio.drops==TRUE) cgf.dt <- cgf.dt[!(abs(age_ratio-1)>drop.criteria)]

#sex ratio: ratio of m:f
#this ratio should be around 1 if the sample is of good quality
cgf.dt[, sex_ratio := sum(sex==1)/sum(sex==2), by=.(nid, source, admin_1)]

#digit preference scores for age, weight, height
#dps is calculated using SMART criteria, with a package called nipnTK
#dps cutoffs based on SMART: https://reliefweb.int/sites/reliefweb.int/files/resources/smart_survey_eastern_ghouta_november_2017_final_.pdf
#Digit preference score: 0-7 excellent, 8-12 good, 13-20 acceptable and > 20 problematic
drop.criteria <- 20

#age
cgf.dt[, age_dps := digitPreference(age_mo, digits=0)$dps, by=.(admin_1, nid)]
message("there are ",
        cgf.dt[age_dps>drop.criteria] %>% uniqueN(by=c('nid', 'source', 'admin_1')),
        " groups that have unreasonable age DPS \n see admin 1s below:")
cgf.dt[age_dps>drop.criteria, paste0(ihme_loc_id, "->", admin_1) %>% table]

#height
cgf.dt[, height_dps := digitPreference(child_height, digits=1)$dps, by=.(admin_1, nid)]
message("there are ",
        cgf.dt[height_dps>drop.criteria] %>% uniqueN(by=c('nid', 'source', 'admin_1')),
        " groups that have unreasonable height DPS \n see admin 1s below:")
#cgf.dt[height_dps>drop.criteria, paste0(ihme_loc_id, "->", admin_1) %>% table]

#weight
cgf.dt[, weight_dps := digitPreference(child_weight, digits=1)$dps, by=.(admin_1, nid)]
message("there are ",
        cgf.dt[weight_dps>drop.criteria] %>% uniqueN(by=c('nid', 'source', 'admin_1')),
        " groups that have unreasonable height DPS \n see admin 1s below:")
#cgf.dt[weight_dps>drop.criteria, paste0(ihme_loc_id, "->", admin_1) %>% table]

#average the 3 indicators
cgf.dt[, mean_dps := apply(.SD, 1, mean), .SDcols=names(cgf.dt)[names(cgf.dt) %like% 'eight_dps'], by=.(admin_1, nid)]
message("there are ",
        cgf.dt[mean_dps>drop.criteria] %>% uniqueN(by=c('nid', 'source', 'admin_1')),
        " groups that have unreasonable average (age/weight/height) DPS \n see admin 1s below:")
#cgf.dt[mean_dps>drop.criteria, paste0(ihme_loc_id, "->", admin_1) %>% table]

#remove the points
if(age.dps.drops==TRUE) cgf.dt <- cgf.dt[!(age_dps>drop.criteria)]
if(mean.dps.drops==TRUE) cgf.dt <- cgf.dt[!(mean_dps>drop.criteria+15)]

#how many groups do we have left?
message('after cgf data cleaning, we are left with this many points:')
uniqueN(cgf.dt, by=c('nid', 'source', 'admin_1'))/og.group.count

#create indicators of various states
cgf.dt[, (indicators) := 0] #begin by filling them in as naught
cgf.dt[, status := "i_none"]
cgf.dt[HAZ > z.cutoff & WAZ > z.cutoff & WHZ > z.cutoff, c('i_none', 'status') := list(1, 'i_none')] #best case
cgf.dt[HAZ > z.cutoff & WAZ > z.cutoff & WHZ <= z.cutoff, c('i_wasting', 'status') := list(1, 'i_wasting')] #just wasted
cgf.dt[HAZ > z.cutoff & WAZ <= z.cutoff & WHZ > z.cutoff, c('i_underweight', 'status')  := list(1, 'i_underweight')] #just underweight
cgf.dt[HAZ <= z.cutoff & WAZ > z.cutoff & WHZ > z.cutoff, c('i_stunting', 'status')  := list(1, 'i_stunting')] #just stunted
cgf.dt[HAZ > z.cutoff & WAZ <= z.cutoff & WHZ <= z.cutoff, c('i_wasting_underweight', 'status')  := list(1, 'i_wasting_underweight')] #both underweight and wasted
cgf.dt[HAZ <= z.cutoff & WAZ <= z.cutoff & WHZ > z.cutoff, c('i_stunting_underweight', 'status')  := list(1, 'i_stunting_underweight')] #both stunted and underweight
cgf.dt[HAZ <= z.cutoff & WAZ <= z.cutoff & WHZ <= z.cutoff, c('i_all', 'status')  := list(1, 'i_all')] #worst case

#also create an indicator for stunted/wasted but not underweight (this is impossible, so dont put it in your vars to loop over)
#this is always 0 but bobby wants it
cgf.dt[, i_stunting_wasting := 0] #begin by filling them in as naught
cgf.dt[HAZ <= z.cutoff & WAZ > z.cutoff & WHZ <= z.cutoff, c('i_stunting_wasting', 'status') := list(1, 'i_stunting_wasting')] #both stunted and wasted (inconceivable!)

#calculate prevalence of each state by admin1
# get weights first
cgf.dt[, pw := weight/sum(weight), by = .(nid, source, admin_1)]
cgf.dt[, (prevs) := lapply(.SD, weighted.mean, w=pw), .SDcols=indicators, by=.(nid, source, admin_1)]
cgf.dt[, p_stunting_wasting := lapply(.SD, weighted.mean, w=pw), .SDcols='i_stunting_wasting', by=.(nid, source, admin_1)]

#summary stats for HAZ, WAZ, WHZ (by admin1)
#mean
cgf.dt[, m_haz := weighted.mean(HAZ, pw), by=.(nid, source, admin_1)]
cgf.dt[, m_waz := weighted.mean(WAZ, pw), by=.(nid, source, admin_1)]
cgf.dt[, m_whz := weighted.mean(WHZ, pw), by=.(nid, source, admin_1)]
cgf.dt[, m_age := weighted.mean(age_mo, pw), by=.(nid, source, admin_1)]
#sd
cgf.dt[, sd_haz := wt.sd(HAZ, pw), by=.(nid, source, admin_1)]
cgf.dt[, sd_waz := wt.sd(WAZ, pw), by=.(nid, source, admin_1)]
cgf.dt[, sd_whz := wt.sd(WHZ, pw), by=.(nid, source, admin_1)]
cgf.dt[, sd_age := wt.sd(age_mo, pw), by=.(nid, source, admin_1)]
#covariance
cgf.dt[, (corr.vars) := cor(.SD) %>% as.list, by=.(nid, source, admin_1), .SDcols=c('HAZ', 'WAZ', 'WHZ')]

#sd for haz
drop.criteria <- 2
message("there are ",
        cgf.dt[sd_haz>drop.criteria] %>% uniqueN(by=c('nid', 'source', 'admin_1')),
        " groups that have a HAZ with SD > 2 \n see admin 1s below:")
cgf.dt[sd_haz>drop.criteria, paste0(ihme_loc_id, "->", admin_1) %>% table]
#drop the points
if(sd.drops==TRUE) cgf.dt<- cgf.dt[sd_haz<=drop.criteria]

message('after cgf data cleaning, we are left with this many points:')
uniqueN(cgf.dt, by=c('nid', 'source', 'admin_1'))/og.group.count


##merge the two indicator data.tables based on the admin1 within each survey
dt <- merge(u5m.dt,
            unique(cgf.dt, by=c('nid', 'source', 'admin_1')),
            by=c('nid', 'source', 'admin_1'))

#do some summarizing and exploration
if (summary.plot==T) {


  #also add in the admin1 level u5m values to cgf dataset for plotting
  cgf.dt <- merge(u5m.dt,
                  cgf.dt,
                  by=c('nid', 'source', 'admin_1'))


  ##xplore data distributions
  #height by status
  ggplot(cgf.dt[year>1999], aes(x=child_height, fill=region_name)) +
    geom_density(alpha=.5) +
    facet_wrap(~status) +
    scale_fill_manual(values=region.palette) +
    theme_minimal()

  #weight by status
  ggplot(cgf.dt[year>1999], aes(x=child_weight, fill=region_name)) +
    geom_density(alpha=.5) +
    facet_wrap(~status) +
    scale_fill_manual(values=region.palette) +
    theme_minimal()

  #weight and height density
  ggplot(cgf.dt[year>1999], aes(x = child_weight, y = child_height, color=status)) +
    geom_point(size=.5, alpha=.1) +
    geom_density_2d() +
    #stat_density2d(aes(fill = ..density..), contour = F, geom = 'tile') +
    facet_wrap(~region_name) +
    scale_color_manual(values=region.palette) +
    theme_minimal()

  #age by u5m quartile
  ggplot(cgf.dt[year>1999], aes(x=age_mo, fill=region_name)) +
    geom_histogram(alpha=.5) +
    facet_wrap(~u5m_bin) +
    scale_fill_manual(values=region.palette) +
    theme_minimal()

}

#logit the respone variable for a guassian regression
dt[, u5m_logit := qlogis(u5m)]
dt <- dt[!dt[is.infinite(u5m_logit)]] #for now we will drop the 20 rows that are 0 for u5m

#normalize the population for u5m to generate model weight
dt[, model_weight := childyrs/mean(childyrs)] #TODO look into this assumption, might not be right to normalize here

#save prepped model data
saveRDS(dt, file=file.path(data.dir, 'prepped', 'model_data_z3.RDS'))

#save a lite version of model data for bobby
lite.vars <- c('super_region_name', 'region_name', 'ihme_loc_id', 'year', 'nid', 'admin_1',
               'u5m', 'childyrs', 'm_haz', 'm_waz', 'm_whz', corr.vars, prevs, 'p_stunting_wasting', 'pop', 'sdi', 'haqi')

saveRDS(dt[, lite.vars, with=F],
        file=file.path(data.dir, 'prepped', 'model_data_lite_z3.RDS'))
write.csv(dt[, lite.vars, with=F],
        file=file.path(data.dir, 'prepped', 'model_data_lite.csv'))
#***********************************************************************************************************************

# ---SUMMARIZE----------------------------------------------------------------------------------------------------------
##calculate coverage and summary statistics for CGF indicators##
if (summary.plot==T) {

##figure 1 scatter##
#new dt
fig1.dt <- copy(cgf.dt)

#save so bobby can play with this figure
#save prepped model data
saveRDS(fig1.dt, file=file.path(data.dir, 'prepped', 'fig1_data.RDS'))

#merge on the admin1 level u5m results
fig1.dt <- merge(fig1.dt,
                 u5m.dt,
                 by=c('nid', 'source', 'admin_1'))

#make a new age/sex variable for the plots
fig1.dt[, age_round := plyr:::round_any(age_mo, 10)]
fig1.dt[sex==1, gender := "Male"]
fig1.dt[sex==2, gender := "Female"]

#create an RGB venn diagram color palette using additive mixing
none.color <- "#B3B3B3"
wasting.color <- "#2B4E76"
underweight.color <- "#F0C814"
stunting.color <-"#DC2828"
wasting.underweight.color <- "#466D37"
stunting.underweight.color <- "#E5681D"
stunting.wasting.color <- "#5D3156"
all.color <- "#83593B"
venn.palette <- c("i_none" = none.color, "i_wasting" = wasting.color, 'i_underweight' = underweight.color, "i_stunting" = stunting.color,
                  'i_wasting_underweight' = wasting.underweight.color, 'i_stunting_underweight' = stunting.underweight.color,
                  'i_stunting_wasting' = stunting.wasting.color, 'i_all'= all.color)

#make the scatterplot and venn diagram legend
plot <- ggplot(data=cgf.dt, aes(x=child_weight, y=child_height, color=age_mo, shape=as.factor(sex))) +
  geom_point() +
  facet_wrap(~status) +
  scale_color_distiller("Age in Months", palette = "YlGnBu") +
  ggtitle("Weight vs Height, by Category") +
  theme_minimal()

print(plot)

ggsave(filename=file.path(graph.dir, 'height_vs_weight_by_category.png'), height=8, width=12, units='in')

plot <- ggplot(data=fig1.dt, aes(x=child_weight, y=child_height, color=status)) +
  geom_point(alpha=.2) +
  facet_wrap(~age_round, nrow=2, ncol=4) +
  scale_color_manual(guide=F, values=venn.palette) +
  scale_alpha_continuous(guide=F) +
  scale_x_continuous("Weight (hg)") +
  scale_y_continuous("Height (mm)") +
  ggtitle("Child Height vs Weight", subtitle="Age groups rounded to the nearest 10 months") +
  theme_minimal()

print(plot)

ggsave(filename=file.path(graph.dir, 'height_vs_weight_by_age.png'), height=8, width=12, units='in')

tiff(file=file.path(graph.dir, 'heigh_vs_weight_by_age_scale.tiff'),
     width=700, height=480)

grid.newpage()
scale <-
  draw.triple.venn(area1 = nrow(fig1.dt[status %like% "stunting|all"]),
                   area2 = nrow(fig1.dt[status %like% "wasting|all"]),
                   area3 = nrow(fig1.dt[status %like% "underweight|all"]),
                   n12 = nrow(fig1.dt[status %like% "stunting_wasting|all"]),
                   n23 = nrow(fig1.dt[status %like% "wasting_underweight|all"]),
                   n13 = nrow(fig1.dt[status %like% "stunting_underweight|all"]),
                   n123 = nrow(fig1.dt[status %like% "all"]),
                   category = c("Stunting", "Wasting", "Underweight"), lty = "blank",
                   fill = c(stunting.color, wasting.color, underweight.color),
                   alpha = c(.6, .5, .5),
                   cex = rep(1.7, 7),
                   cat.cex = rep(2.5, 3),
                   cat.pos = c(330, 30, 0),
                   cat.dist = c(.075, .075, .02),
                   print.mode = "raw")

dev.off()


##mapping##
mapPoints(cgf.dt, 2010, 2016, 'p_all', -1)

##scatters##
scatterVars(dt, 2000, 2016, 'p_all', 'u5m', 'region_name')

pdf(file=file.path(graph.dir, 'u5m_vs_prevalence_scatter_2000_2016.pdf'),
                   width=8, height=8)

lapply(prevs, scatterVars,
       plot.dt=dt, year.start=2000, year.end=2016, y.var='u5m', color.var='region_name') %>%
  print

dev.off()

pdf(file=file.path(graph.dir, 'u5m_vs_prevalence_scatter_1990_2000.pdf'),
    width=8, height=8)

lapply(prevs, scatterVars,
       plot.dt=dt, year.start=1990, year.end=2000, y.var='u5m', color.var='region_name') %>%
  print

dev.off()

}
#***********************************************************************************************************************

# ---MODELLING----------------------------------------------------------------------------------------------------------

##graph predictions##
pdf(file=file.path(graph.dir, 'z3_predictions_2000_2016.pdf'),
    width=8, height=8)

##create some GLMs and GAMs of u5m against the indicator prevs##
preds.guass <- mclapply(prevs, gamPreds,
                        input.dt=dt, model.fam='gaussian',
                        year.start=2000, year.end=2016,
                        response.var='u5m_logit',
                        k.value=4, by.var='region_name', weight.var='childyrs',
                        mc.cores=cores) %>%
  rbindlist

scatterPreds(plot.dt=preds.guass, year.start=2000, year.end=2016, color.var='region_name',
             this.title="U5M (probability) Predictions, ",
             this.subtitle="Guassian family, with Region Smoothers & Population-Weighting")

##try a super region model

preds.guass <- mclapply(prevs, gamPreds,
                        input.dt=dt, model.fam='gaussian',
                        year.start=2000, year.end=2016,
                        response.var='u5m_logit',
                        k.value=4, by.var='super_region_name', weight.var='childyrs',
                        mc.cores=cores) %>%
  rbindlist

scatterPreds(plot.dt=preds.guass, year.start=2000, year.end=2016, color.var='super_region_name',
             this.title="U5M (probability) Predictions, ",
             this.subtitle="Guassian family, with Super Region Smoothers & Population-Weighting")

##try a global model
preds.guass <- mclapply(prevs, gamPreds,
                        input.dt=dt, model.fam='gaussian',
                        year.start=2000, year.end=2016,
                        response.var='u5m_logit',
                        k.value=4, by.var='source', weight.var='childyrs',
                        mc.cores=cores) %>%
  rbindlist


scatterPreds(plot.dt=preds.guass, year.start=2000, year.end=2016, color.var='source',
             this.title="U5M (probability) Predictions, ",
             this.subtitle="Guassian family, with Global Smoothers & Population-Weighting")


dev.off()

#***********************************************************************************************************************

# ---SCRAP--------------------------------------------------------------------------------------------------------------
if(run.scrap==T){
preds.wt.quasi <- mclapply(prevs, gamPreds,
                           input.dt=dt, model.fam='quasibinomial',
                           year.start=2000, year.end=2016,
                           response.var='u5m',
                           k.value=4, by.var='region_name', weight.var='model_weight',
                           mc.cores=cores) %>%
  rbindlist

preds.wt <- mclapply(prevs, gamPreds,
                     input.dt=dt, model.fam='binomial',
                     year.start=2000, year.end=2016,
                     response.var='u5m',
                     k.value=4, by.var='region_name', weight.var='model_weight',
                     mc.cores=cores) %>%
  rbindlist

preds.wt.betar <- mclapply(prevs, gamPreds,
                           input.dt=dt, model.fam='betar',
                           year.start=2000, year.end=2016,
                           response.var='u5m',
                           k.value=4, by.var='region_name', weight.var='model_weight',
                           mc.cores=cores) %>%
  rbindlist

scatterPreds(plot.dt=preds.wt.quasi, year.start=2000, year.end=2016, color.var='region_name',
             this.title="U5M (probability) Predictions, ",
             this.subtitle="Quasibinomial family, with Region Smoothers & Population-Weighting")

scatterPreds(plot.dt=preds.wt, year.start=2000, year.end=2016, color.var='region_name',
             this.title="U5M (probability) Predictions, ",
             this.subtitle="Binomial family, with Region Smoothers & Population-Weighting")

scatterPreds(plot.dt=preds.wt.betar, year.start=2000, year.end=2016, color.var='region_name',
             this.title="U5M (probability) Predictions, ",
             this.subtitle="BetaR family, with Region Smoothers & Population-Weighting")
}
