#
# Build the pollen time series data file using default age models in Neotoma or Calcote data set
#

library(neotoma)

source('r/utils/make_data_funs.r')
source('r/utils/compile_lists.r')

# translation table to convert from aggregated taxa to stepps taxa
pollen.equiv.stepps = read.csv("pollen.equiv.stepps.csv", stringsAsFactors=F)

new_meta = FALSE

###############################################################################################################

if (new_meta){  
  # get new metadata
  pollen_meta_neo <- get_neo()
  pollen_meta_neo <- pollen_meta_neo[with(pollen_meta_neo, order(handle)),]
  pollen_meta_neo <- split_mi(pollen_meta_neo, longlat=TRUE)
  write.table(pollen_meta_neo, paste('data/pollen_meta_', Sys.Date(), '.csv', sep=''), row.names=FALSE, sep=',', na='')
} else {
  # for the paper based on the elicitation exercise
  pollen_meta_neo <- read.csv('data/pollen_meta_2014-07-22.csv', header=TRUE, stringsAsFactors=FALSE, sep=',')
}

nsites = nrow(pollen_meta_neo)

# download and save the raw data
pollen2k_raw = list()
for (i in 1:nsites){ 
  print(i)
  
  id = pollen_meta_neo$datasetID[i]
  if (id == 1394) id = 15274
  
  pollen2k_raw[[i]] = get_download(id)
}
save(pollen2k_raw, file=paste0('data/pollen2k_raw_', Sys.Date(), '.rdata'))

# aggregate the taxa
pollen2k = list()
for (i in 1:nsites){  
  print(i)
  convert1 = compile_list_neotoma(pollen2k_raw[[i]], 'Stepps')
  #pollen2k[[i]] = compile_list_stepps(convert1, list.name='all', pollen.equiv.stepps, cf = TRUE, type = TRUE)
  pollen2k[[i]] = compile_list_stepps(convert1, list.name='must_have', pollen.equiv.stepps, cf = TRUE, type = TRUE)
  pollen2k[[i]]$dataset$site.data$state = pollen_meta_neo$state[i]
}

save(pollen2k, file=paste0('data/pollen2k_', Sys.Date(), '.rdata'))

###############################################################################################################
# Build the ts data
###############################################################################################################

pollen_ts = list()
age_models = vector(length=0)
types = vector(length=0)

for (i in 1:nsites){
  
  print(i)
  
  id       <- pollen2k[[i]]$dataset$dataset.meta$dataset.id
  lat      <- pollen2k[[i]]$dataset$site.data$lat
  long     <- pollen2k[[i]]$dataset$site.data$long
  altitude <- pollen2k[[i]]$dataset$site.data$elev
  state    <- pollen2k[[i]]$dataset$site.data$state
  handle   = x$dataset$dataset.meta$collection.handle
  sitename <- gsub("[ ']", "", as.vector(pollen2k[[i]]$dataset$site.data$site.name))
  
  chrons_full   = sapply(pollen2k[[i]]$sample.meta[1, 'chronology.name'], as.character)
  
  chrons_split = strsplit(chrons_full, split=" ")
  chron        = unname(sapply(chrons_split, "[[", 1))
  type         = sapply(pollen2k[[i]]$sample.meta[1, 'age.type'], as.character)
  ages         = pollen2k[[i]]$sample.meta$age
  
  age_models = c(age_models, chron) 
  types      = c(types, type)
  print(chron)
  
  meta = list(id = id, sitename = sitename, lat = lat, long = long, state = state, altitude = altitude, model = chron)
  
  pollen_ts[[i]] = list(meta = meta, counts = data.frame(ages, pollen2k[[i]]$counts))
  #   pollen_ts_nap[[i]] = list(meta = meta, counts = data.frame(ages, pollen2k[[i]]$counts_pls, pollen2k[[i]]$counts_nap))
}


###############################################################################################################
# Add Calcote cores that have not been added to Neotoma
###############################################################################################################

long_cores <- c('Ferr01X', 'Hell Hole', 'Lone02', 'Warner')
handles    <- c('FERRY', 'HELLHOLE', 'LONE', 'WARNER')
ids        <- c('CLH2', 'CLH3', 'CLH5', 'CLH6')

clh_meta   <- read.csv('data/hotchkiss_lynch_calcote_meta_v0.1.csv', header=TRUE)
clh_counts <- read.csv('data/hotchkiss_lynch_calcote_counts_v0.1.csv', header=TRUE, stringsAsFactors=FALSE)

clh_meta <- clh_meta[clh_meta$name %in% long_cores, ]
clh_counts <- clh_counts[clh_counts$name %in% long_cores, ]

for (i in 1:length(long_cores)){
  
  print(i)
  
  id       <- ids[i]
  lat      <- clh_meta$lat[i]
  long     <- clh_meta$long[i]
  altitude <- clh_meta$Elev.m.[i]
  state    <- 'wisconsin'
  sitename <- long_cores[i]
  chron    <- 'clh'
  
  meta = list(id = id, sitename = sitename, lat = lat, long = long, state = state, altitude = altitude, 
              model = chron)
  
  ages   = clh_counts[which(clh_counts$name == sitename), 'age']
  counts = clh_counts[which(clh_counts$name == sitename), 8:ncol(clh_counts)]
  
  convert1      = compile_list_neotoma(counts, 'Stepps')
  counts_stepps = compile_list_stepps(convert1, list.name='must_have', pollen.equiv.stepps, cf = TRUE, type = TRUE)
  
  pollen_ts[[nsites+i]] = list(meta = meta, counts = data.frame(ages, counts_stepps))
  #   pollen_ts_nap[[i]] = list(meta = meta, counts = data.frame(ages, pollen2k[[i]]$counts_pls, pollen2k[[i]]$counts_nap))
}


###############################################################################################################
# Add Calcote cores that have not been added to Neotoma
###############################################################################################################

#write it to a csv file: pollen_ts.csv
pollen_ts_df = list()
for (i in 1:length(pollen_ts)){
  print(i)
  counts = pollen_ts[[i]]$counts
  meta   = t(replicate(nrow(counts),unlist(pollen_ts[[i]]$meta)))
  pollen_ts_df = rbind(pollen_ts_df, cbind(meta, counts))
}
write.table(pollen_ts_df, file=paste('data/pollen_ts_', Sys.Date(), '.csv', sep=''), quote=FALSE, row.names=FALSE)
