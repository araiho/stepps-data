# build a calibration data set based on our pre-settlement exercise results
# cal.min: uses lowest depth assigned by an expert as pre-settlement
# cal.max: uses highest depth assigned by an expert as pre-settlement

#####################################################################################
# user pars
#####################################################################################
library(neotoma)

# source('~/Documents/paleon/stepps2/eda/r/utils/compile_list_stepps.r')
source('r/utils/compile_lists_v2.r')
source('r/utils/build_data_funs.r')

# translation table to convert from aggregated taxa to stepps taxa
pollen.equiv.stepps = read.csv("pollen.equiv.stepps.csv", stringsAsFactors=F)
pollen.equiv = read.csv("pollen.equiv.csv", stringsAsFactors=F)

elicit = read.csv('data/pollen_meta_2014-05-01_compiled.csv', header=TRUE, stringsAsFactor=FALSE)
elicit$datasetID[which(elicit$datasetID == 1394)] = 15274

version = 'v3'

list_name = 'must_have'
# list_name = 'kujawa'

#####################################################################################
# metadata for calibration paper
#####################################################################################

meta = data.frame(Site=character(0), ID=character(0),  Lat=numeric(0), Long=numeric(0),PI=character(0), 
                  Depth=numeric(0), Neotoma=character(0), Calibration=character(0), Notes=character(0))

#####################################################################################
# read in dictionary and pollen data
#####################################################################################
neo_ids = as.numeric(elicit$datasetID[!(substr(elicit$datasetID, 1, 3) == 'CLH')])
clh_ids = elicit$datasetID[substr(elicit$datasetID,1,3) == 'CLH' ]

# neo_ids[which(neo_ids == 1394)] = 15274

# load list containing pollen counts
# will load object called pollen_dat
# the first time takes a while to pull from Neotoma
if (!file.exists('data/pollen_neo.rdata')){
  pollen_neo <- get_download(neo_ids)
  save(pollen_neo, file='data/pollen_neo.rdata')
} else {
  load('data/pollen_neo.rdata') 
}

taxa = sort(unique(pollen.equiv.stepps$must_have))

# split ids
# get rid of mi:south for now
# elicit = elicit[elicit$state != 'michigan:south', ]

#####################################################################################
# user pars
#####################################################################################

# number of site where all were in agreement
agree = apply(elicit[,9:12], 1, function(x) length(unique(x)))
sum(agree==1)
sum(agree==2)
sum(agree==3)
sum(agree==4)

elicit[agree==1, c(1:3, 9:16)]

elicit[which(agree==4), 'handle']
elicit[agree==4, c(1:3, 9:16)]

no_pre = apply(elicit[,9:12], 1, function(x) all(x==-1))
sum(no_pre)
elicit[which(no_pre),1:3]

elicit_out = elicit[,c(1,9:16)]
elicit_out[is.na(elicit_out)] = 0

#####################################################################################
# build neotoma calibration data sets from elicitation results
#####################################################################################
#(data, ids, elicit, taxa, type, depth_type)
data=pollen_neo
ids=neo_ids
elicit=elicit
type='neo'
depth_type='mid'
meta=meta

# build the neo calibration datasets

# function(data, ids, meta, taxa, type, depth_type, meta){
# neo_cal_min_out = get_cal_data(pollen_neo, neo_ids, elicit, taxa, type='neo', depth_type='min', meta)
# neo_cal_max_out = get_cal_data(pollen_neo, neo_ids, elicit, taxa, type='neo', depth_type='max', meta)
neo_cal_mid_out = get_cal_data(pollen_neo, neo_ids, elicit, type='neo', depth_type='mid', meta, list_name)
# neo_cal_min = neo_cal_min_out$cal
# neo_cal_max = neo_cal_max_out$cal
neo_cal_mid = neo_cal_mid_out$cal

# save the rdata
# save(neo_cal_min, file=paste0('data/cal_min_', version, '.rdata'))
# save(neo_cal_max, file=paste0('data/cal_max_', version, '.rdata'))
# save(neo_cal_mid, file=paste0('data/cal_mid_', version, '.rdata'))

# function(data, ids, meta, taxa, type, depth_type, meta){
# neo_cal_min_out = get_cal_data(pollen_neo, neo_ids, elicit, taxa, type='neo', depth_type='min', meta)
# neo_cal_max_out = get_cal_data(pollen_neo, neo_ids, elicit, taxa, type='neo', depth_type='max', meta)
neo_cal_mid_out = get_cal_data(pollen_neo, neo_ids, elicit, type='neo', depth_type='mid', meta, list_name)
# neo_cal_min = neo_cal_min_out$cal
# neo_cal_max = neo_cal_max_out$cal
neo_cal_mid = neo_cal_mid_out$cal

meta_neo = neo_cal_mid_out$meta

#####################################################################################
# build calcote calibration data sets from elicitation results

clh <- read.csv('data/pollen_counts_calcote.csv', stringsAsFactors=FALSE)

data=clh
ids=clh_ids
elicit=elicit
type='clh'
depth_type='mid'
# 
# # function(data, ids, meta, taxa, type, depth_type)
# clh_cal_min_out = get_cal_data(clh, clh_ids, elicit, taxa, type='clh', depth_type='min', meta, list_name)
# clh_cal_max_out = get_cal_data(clh, clh_ids, elicit, taxa, type='clh', depth_type='max', meta, list_name)
# clh_cal_mid_out = get_cal_data(clh, clh_ids, elicit, taxa, type='clh', depth_type='mid', meta, list_name)
# clh_cal_min = clh_cal_min_out$cal
# clh_cal_max = clh_cal_max_out$cal
# clh_cal_mid = clh_cal_mid_out$cal


# function(data, ids, meta, taxa, type, depth_type)
clh_cal_min_out = get_cal_data(clh, clh_ids, elicit, type='clh', depth_type='min', meta, list_name)
clh_cal_max_out = get_cal_data(clh, clh_ids, elicit, type='clh', depth_type='max', meta, list_name)
clh_cal_mid_out = get_cal_data(clh, clh_ids, elicit, type='clh', depth_type='mid', meta, list_name)
clh_cal_min = clh_cal_min_out$cal
clh_cal_max = clh_cal_max_out$cal
clh_cal_mid = clh_cal_mid_out$cal

meta_clh = clh_cal_mid_out$meta

#####################################################################################
# build calibration data set from sites with only core-top and pre-settlement
#####################################################################################

clh.sites  <- read.csv('data/hotchkiss_lynch_calcote_meta_v0.1.csv', stringsAsFactors = FALSE, sep=',')
clh.counts <- read.csv('data/hotchkiss_lynch_calcote_counts_v0.1.csv', stringsAsFactors = FALSE, sep=',')

clh.sites$name <- gsub(" ","", clh.sites$name, fixed=TRUE)
clh.counts$name <- gsub(" ","", clh.counts$name, fixed=TRUE)

get_clhpre <- function(clh.sites, clh.counts, meta){
  
  n <- nrow(clh.sites)
  
  ids    <- vector(mode="numeric", length=0)
  states <- vector(mode="numeric", length=0)
  lat    <- vector(mode="numeric", length=0)
  long   <- vector(mode="numeric", length=0)
  site   <- vector(mode="character", length=0)
  pi     <- vector(mode="character", length=0)
  pre.all    <- matrix(NA, nrow=0, ncol=length(taxa))
  
  site.count = 0
  
  for (i in 1:n){
    
    site.i <- as.character(clh.sites$name[i])
    idx    <- which(clh.counts$name == site.i)
    sample.type   <- clh.counts[idx,1]
    counts = clh.counts[idx,]
    
    if ((length(sample.type) <= 2) & (any(sample.type=='Pre'))){
      
      pre.site = as.matrix(counts[which(sample.type=="Pre"),8:ncol(clh.counts)])
      compressed2stepps = compile_list_neotoma(pre.site, 'Stepps')
      compressed2model = compile_list_stepps(compressed2stepps, list.name=list_name, pollen.equiv.stepps, cf = TRUE, type = TRUE)
      #       compressed2model  = compile_taxa_stepps(compressed2stepps, list.name='must_have', 
      #                                               alt.table=pollen.equiv.stepps, cf = TRUE, type = TRUE)
      #       compressed2model  = t(compressed2model[,!(colnames(compressed2model)=='Other')])
      pre.all = rbind(pre.all, compressed2model)
      site.count = site.count + 1
      ids    = c(ids, paste('CALPRE', site.count, sep=''))
      states = c(states, 'wisconsin')
      site   = c(site, site.i)
      lat    = c(lat, clh.sites$lat[i])
      long   = c(long, clh.sites$long[i])
      pi     = c(pi, clh.counts$analyst[idx[1]])
      
      meta.row = data.frame(Site     = site.i, 
                            ID       = paste0('CLHPRE', site.count),
                            Lat      = clh.sites$lat[i], 
                            Long     = clh.sites$long[i], 
                            PI       = clh.counts$analyst[idx[1]], 
                            Depth    = counts[sample.type=='Pre', 'depth_mid'],
                            Neotoma  = 'N', 
                            Calibration = 'Y',
                            Notes    = '')
      meta  = rbind(meta, meta.row)
    } else {
      next
    }
    
  }
  
  clh_pre <- data.frame(id=ids, site=site, long=long, lat=lat, state=states, pi=pi, depth=rep(NA, length(ids)), pre.all, 
                        stringsAsFactors=FALSE)
  
  return(list(clh_pre=clh_pre, meta_clhpre = meta))
  
}

clh_pre_out = get_clhpre(clh.sites, clh.counts, meta)
clh_pre     = clh_pre_out$clh_pre
meta_clhpre = clh_pre_out$meta_clhpre

#########################################################################################################################################
# build calibration data set from sites with only core-top and pre-settlement
#########################################################################################################################################

# cal_min = rbind(neo_cal_min, clh_cal_min, clh_pre)
# cal_max = rbind(neo_cal_max, clh_cal_max, clh_pre)
cal_mid = rbind(neo_cal_mid, clh_cal_mid, clh_pre)

# save tha data files
suffix = version
if (list_name == 'kujawa'){
  suffix = paste0(list_name, '_', suffix)
  colnames(cal_mid)[which(colnames(cal_mid) == 'CEDAR.JUNIPER')] = 'JUNIPER'
  colnames(cal_mid)[which(colnames(cal_mid) == 'POPLAR.TULIP.POPLAR')] = 'POPLAR'
  colnames(cal_mid) = tolower(colnames(cal_mid))
}

# write.table(cal_min, file=paste('data/cal_data_upper_depth_', Sys.Date(), '.csv', sep=''), sep=',', row.names=FALSE)
# write.table(cal_max, file=paste('data/cal_data_lower_depth_', Sys.Date(), '.csv', sep=''), sep=',', row.names=FALSE)
write.table(cal_mid, file=paste('data/cal_data_mid_depth_', suffix, '.csv', sep=''), sep=',', row.names=FALSE)

meta      = rbind(meta_neo, meta_clh, meta_clhpre)
col_names = paste0('{', colnames(meta), '}')

# csv to latex package doesn't like NA
meta[is.na(meta)] = ''
# use { and } to indicate table entries
meta = t(apply(meta, 1, function(x) paste0('{', x, '}')))
colnames(meta) = col_names

write.table(meta, file=paste('data/meta_', version, '.csv', sep=''), sep=',', row.names=FALSE, quote=FALSE)
