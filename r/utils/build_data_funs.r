get_cal_data <- function(data, ids, elicit, taxa, type, depth_type, meta){
  
  cal = data.frame(matrix(NA, nrow=0, ncol=length(taxa)+7))
  # colnames(cal) = taxa
  
  ages_diff = data.frame(matrix(NA, nrow=0, ncol=4))
  
  skipped.id.3na=list()
  skipped.site.3na=list()
  
  skipped.id.2na=list()
  skipped.site.2na=list()
  
  if (type == 'neo'){
    ids.data = vector(length=length(data))
    for (i in 1:length(data)){
      ids.data[i] = data[[i]]$dataset$dataset.meta$dataset.id
    }
  }
  
  for (i in 1:length(ids)){  
    print(i)
    
    if (type=='neo'){
      id      = as.numeric(ids[i])
      from_db = 'Y' 
      
      print(id)
      
      # skip green lake
      if (id == 982){  
        print('Skipping Green Lake')
        meta.row = data.frame(site     = site, 
                              id       = id, 
                              long     = long, 
                              lat      = lat, 
                              pi       = pi, 
                              depth    = depth, 
                              neotoma  = from_db, 
                              calibration = 'N',
                              notes    = 'C')
        meta  = rbind(meta, meta.row)
        next 
      }
      
      # skip ryerse lake
      if (id == 2318){  
        print('Skipping Ryerse Lake')
        meta.row = data.frame(site     = site, 
                              id       = id, 
                              long     = long, 
                              lat      = lat, 
                              pi       = pi, 
                              depth    = depth, 
                              neotoma  = from_db, 
                              calibration = 'N',
                              notes    = 'C')
        meta  = rbind(meta, meta.row)
        next 
      }
      
      # skip lake mary
      if (id == 3473){ 
        print('Skipping Lake Mary')
        meta.row = data.frame(site     = site, 
                              id       = id, 
                              long     = long, 
                              lat      = lat, 
                              pi       = pi, 
                              depth    = depth, 
                              neotoma  = from_db, 
                              calibration = 'N',
                              notes    = 'C')
        meta  = rbind(meta, meta.row)
        next 
      } 
      
      # skip chippewa lake
      if (id == 369){ 
        print('Skipping Chippewa Bog')
        meta.row = data.frame(site     = site, 
                              id       = id, 
                              long     = long, 
                              lat      = lat, 
                              pi       = pi, 
                              depth    = depth, 
                              neotoma  = from_db, 
                              calibration = 'N',
                              notes    = 'C')
        meta  = rbind(meta, meta.row)
        next 
      } 
      
      idx.meta = which(elicit$datasetID == id)      
      id.meta  =  elicit[idx.meta, 'datasetID']
      
      idx.data = which(ids.data == id)
      id.data = data[[idx.data]]$dataset$dataset.meta$dataset.id
      
      site  = elicit[idx.meta, 'site']
      long  = elicit[idx.meta, 'long']
      lat   = elicit[idx.meta, 'lat']
      state = elicit[idx.meta, 'state']
      pi    = elicit[idx.meta, 'pi']
      
      print(site)
      
      if ((id != id.meta) | (id != id.data)){
        stop('Naming or ordering inconsistency between meta file and loaded rdata file!')
      }
      
      x = compile_list_neotoma(data[[idx.data]], 'Stepps')
      x = compile_list_stepps(x, list.name='must_have', pollen.equiv.stepps, cf = TRUE, type = TRUE)
      
      #       x = compile_taxa_stepps(data[[i]], list.name='WhitmoreSmall', alt.table=NULL, cf=TRUE, type=TRUE)
#       x = compile_taxa_stepps(x, list.name='must_have', alt.table=pollen.equiv.stepps, cf = TRUE, type = TRUE)
#       x$counts = x$counts[,!(colnames(x$counts)=='Other')]
#       
#       # add empty columns for missing taxa
#       zero_taxa = taxa[!(taxa %in% colnames(x$counts))]
#       add_back   = matrix(0, nrow=nrow(x$counts), ncol=length(zero_taxa))
#       colnames(add_back) = zero_taxa
#       
#       tmp      = cbind(x$counts, add_back)
#       x$counts = tmp[, sort(colnames(tmp))]
      
      age = x$sample.meta$age
      
      good_samps  = which(age < 2000)
      good_ages   = age[good_samps] 
      good_depths = x$sample.meta$depth[good_samps]
      
      pol = x$counts[good_samps,]
      
    } else if (type=='clh'){
      id = ids[i]
      from_db = 'N'
      
      idx = which(data$datasetID == id)
      
      site  = data$site[idx[1]]
      lat   = data$lat[idx[1]]
      long  = data$long[idx[1]]
      state = data$state[idx[1]]
      pi    = data$pi[idx[1]]
      
      x      = data[idx,1:ncol(data)]
      depths = data$depth_mid[idx]
      
      clh_meta = x[,1:15]
      counts   = x[,16:ncol(x)]
      
      x = compile_list_neotoma(counts, 'Stepps')
      x = compile_list_stepps(x, list.name='must_have', pollen.equiv.stepps, cf = TRUE, type = TRUE)
      
      no_age = any(is.na(clh_meta$age))
      
      if (no_age){
        age = rep(NA, length=nrow(x))
      } else{
        age = clh_meta$age
      }
      
      good_samps  <- which(age < 2000)
      good_ages   <- age[good_samps] 
      good_depths <- depths[good_samps]      
      
      pol <- x[good_samps,]
  
    }
    
    if (!all(taxa == colnames(pol))) {
      print("Taxon mismatch in pollen records.")
      print(paste("--------> site: ", id, sep=''))
    }
    
    # get the results from the exercise
    sh_depths = elicit[elicit$datasetID == id, 9:12]
    sh_depths[sh_depths<=0] = NA
    sh_unsure = sum(elicit[elicit$datasetID == id, 13:16], na.rm=TRUE)
    
    ages_diff_site = rep(NA, 4)
    names(ages_diff_site) = c('pre1', 'pre2', 'pre3', 'pre4')
    for (j in 1:length(sh_depths)){  
      depth = as.numeric(sh_depths[j])
      if (!is.na(sh_depths[j]))
        ages_diff_site[j] = 100 - good_ages[depth]
    }
    
    na_count = sum(is.na(ages_diff_site))
    
    ages_diff_row = data.frame(id = id, site=site, lat=lat, lon=long, num_depths=length(good_ages), 
                              sh_unsure=sh_unsure, na_count=na_count, t(ages_diff_site)) 
    ages_diff     = rbind(ages_diff, ages_diff_row)  
    
    # skip sites for which 3 or 4 experts did not assign a pre sample
    if (na_count >=3 ){
      print(i)
      skipped.id.3na   = c(skipped.id.3na, id)
      skipped.site.3na = c(skipped.site.3na, site)
      meta.row = data.frame(site     = site, 
                            id       = id, 
                            long     = long, 
                            lat      = lat, 
                            pi       = pi, 
                            depth    = depth, 
                            neotoma  = from_db, 
                            calibration = 'N',
                            notes    = 'A')
      meta  = rbind(meta, meta.row)
      next
    } else if ((na_count==2) & (sum(abs(ages_diff_site)>300, na.rm=TRUE)==2)) {
      print(i)
      skipped.id.2na   = c(skipped.id.2na, id)
      skipped.site.2na = c(skipped.site.2na, site)
      meta.row = data.frame(site     = site, 
                            id       = id, 
                            long     = long, 
                            lat      = lat, 
                            pi       = pi, 
                            depth    = depth, 
                            neotoma  = from_db, 
                            calibration = 'N',
                            notes = 'B')
      meta  = rbind(meta, meta.row)
      next
    } else if (id == 1548){ # kotirant
      pre_samp = 1
    } else {
      
      # sort depths
      sh_depths  = sort(sh_depths)
      
      if (depth_type=='min') {
        pre_samp = min(sh_depths, na.rm=TRUE) #sh_depths[which.min(sh_depths)]
      } else if (depth_type=='max') {
        pre_samp = max(sh_depths, na.rm=TRUE) #which.max(sh_depths)
      } else if (depth_type=='mid') {
        pre_samp = as.numeric(sh_depths[length(sh_depths)-1])
      } else {
        stop('Depth type is not one of min, mid, or max.')
      }
    }

    names(pol[pre_samp,]) = colnames(pol)    

    if (type=='neo'){
      pol = t(pol[pre_samp,])
    } else if (type=='clh'){
      #pol = t(pol[pre_samp,])
      pol = pol[pre_samp,]
    }
    
    depth = good_depths[pre_samp]
    
    print(i)
    cal   = rbind(cal, data.frame(id, site, long, lat, state, pi, depth=depth, pol))
    meta.row = data.frame(site     = site, 
                          id       = id, 
                          long     = long, 
                          lat      = lat, 
                          pi       = pi, 
                          depth    = depth, 
                          neotoma  = from_db, 
                          calibration = 'Y',
                          notes    = '')
    meta  = rbind(meta, meta.row)
    
  }
  
  skipped3na = data.frame(id=unlist(skipped.id.3na), site=unlist(skipped.site.3na))  
  skipped2na = data.frame(id=unlist(skipped.id.2na), site=unlist(skipped.site.2na)) 
  
  return(list(cal=cal, meta=meta, skipped3na=skipped3na, skipped2na=skipped2na))
  
}