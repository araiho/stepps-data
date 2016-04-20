data.dir = c("~/stepps-data/data/")
fig.dir = c('~/stepps-data/figures/')

#####
##### Prediction Data #####
#####

##### Install #####
 library(gtools)
 library(neotoma)
# library(reshape)
# library(ggplot2)
# require(grid)
 require(plyr)
 require(maps)
# library(gtools)
# library(sp)
# library(rgdal)
# library(fields)
# library(maptools)
# library(mapplots)
# gpclibPermit()


#####
##### Download Data for MN and WI and MI #####
#####
gpids <- get_table(table.name='GeoPoliticalUnits')
gpid <- gpids[which(gpids$GeoPoliticalName == 'Minnesota'),1]
gpid1 <- gpids[which(gpids$GeoPoliticalName == 'Wisconsin'),1]
gpid2 <- gpids[which(gpids$GeoPoliticalName == 'Michigan'),1]
meta <- get_dataset(datasettype='pollen', gpid=c(gpid), ageyoung=0) #Pollen data for all of MN&WI
meta1 <- get_dataset(datasettype='pollen', gpid=c(gpid1), ageyoung=0) #Pollen data for all of MN&WI
meta2 <- get_dataset(datasettype='pollen', gpid=c(gpid2), ageyoung=0) #Pollen data for all of MN&WI

#loc.PI = rep(0,length(c(meta)))
#for(i in 1:length(loc.PI)) loc.PI[i] = as.character(meta[[i]]$pi.data$ContactName)

hk_counts = read.csv(paste0(data.dir,"hotchkiss_lynch_calcote_counts_v0.1.csv"))
hk_meta = read.csv(paste0(data.dir,"hotchkiss_lynch_calcote_meta_v0.1.csv"))

site.locs <- ldply(meta, function(x) c(x$site.data$long, x$site.data$lat))
site.locs1 <- ldply(meta1, function(x) c(x$site.data$long, x$site.data$lat))
site.locs2 <- ldply(meta2, function(x) c(x$site.data$long, x$site.data$lat))
site.locs<-rbind(site.locs,site.locs1,site.locs2)

pdf(paste0(fig.dir,"all.sites.neotoma.hotchkiss.pdf"))
map('state', xlim=range(site.locs$V1)+c(-2, 2), ylim=range(site.locs$V2)+c(-1, 1))
points(site.locs$V1, site.locs$V2, pch=19, cex=1,col="black")
points(hk_meta$long, hk_meta$lat, pch=19, cex=1,col="blue")
title("All Pollen Sites")
dev.off()

mnwi1 <- rep(0,length(c(meta)))
mnwi2 <- rep(0,length(c(meta1)))
mnwi3 <- rep(0,length(c(meta2)))
for(i in 1:length(meta)) mnwi1[i] <- meta[[i]]$dataset.meta$dataset.id
for(i in 1:length(meta1)) mnwi2[i] <- meta1[[i]]$dataset.meta$dataset.id
for(i in 1:length(meta2)) mnwi3[i] <- meta2[[i]]$dataset.meta$dataset.id
datasets <- as.vector(c(mnwi1,mnwi2,mnwi3))
dat.mnwi <- lapply(datasets, function(x)try(get_download(x)))#get_download(x = as.vector(c(mnwi1,mnwi2,mnwi3)))
#save(dat.mnwi,file=paste0("neotoma_pollen_dat_babystepps",Sys.Date(),".rdata"))
load(paste0(data.dir,file="neotoma_pollen_dat_babystepps2016-04-19.rdata")) # start here unless you think there might be new data in neotoma

##### prepare data frame for adding counts
source("~/stepps-data/r/utils/compile_lists_v2.r")

x=dat.mnwi

# translation table to convert from aggregated taxa to stepps taxa
pollen.equiv.stepps = read.csv("pollen.equiv.stepps.csv", stringsAsFactors=F)
pollen.equiv = read.csv("pollen.equiv.csv", stringsAsFactors=F)

babe.dat <- compile_list_stepps(compile_list_neotoma(dat.mnwi[[1]][[1]],"WhitmoreSmall"),
                    'babystepps', pollen.equiv.stepps = pollen.equiv.stepps)

pol.cal.count<-data.frame(cbind(rep(babe.dat$dataset$site$site.id,length(babe.dat$sample.meta$age)),
rep(babe.dat$dataset$site$lat,length(babe.dat$sample.meta$age)),
    rep(babe.dat$dataset$site$long,length(babe.dat$sample.meta$age)),
        rep(babe.dat$dataset$dataset.meta$dataset.id,length(babe.dat$sample.meta$age)),
            rep(babe.dat$dataset$dataset.meta$collection.handle,length(babe.dat$sample.meta$age)),
                babe.dat$sample.meta$age,babe.dat$counts))

for(i in 2:length(dat.mnwi)){
  babe.dat <- compile_list_stepps(compile_list_neotoma(dat.mnwi[[i]][[1]],"WhitmoreSmall"),
                                  'babystepps', pollen.equiv.stepps = pollen.equiv.stepps)
  pol.cal.count.save<-data.frame(cbind(rep(babe.dat$dataset$site$site.id,length(babe.dat$sample.meta$age)),
                                  rep(babe.dat$dataset$site$lat,length(babe.dat$sample.meta$age)),
                                  rep(babe.dat$dataset$site$long,length(babe.dat$sample.meta$age)),
                                  rep(babe.dat$dataset$dataset.meta$dataset.id,length(babe.dat$sample.meta$age)),
                                  rep(babe.dat$dataset$dataset.meta$collection.handle,length(babe.dat$sample.meta$age)),
                                  babe.dat$sample.meta$age,babe.dat$counts))
  pol.cal.count <- rbind(pol.cal.count,pol.cal.count.save)
}


head(pol.cal.count)

##### Fix problems with matrix
rownames(pol.cal.count)<-seq(1,nrow(pol.cal.count),1)
pol.cal.count[is.na(pol.cal.count)]<-0
colnames(pol.cal.count)<-c("SiteID","LatitudeNorth","LongitudeWest","dataset.id","ContactName","Age",colnames(pol.cal.count[,7:ncol(pol.cal.count)]))
#save(pol.cal.count,file=paste0(data.dir,"pol.cal.count.mnwi3.csv"))

load(paste0(data.dir,"pol.cal.count.mnwi3.csv"))

### CONVERT HOTCHKISS SITES
hk_counts[is.na(hk_counts)] <- 0

hk_counts1 = rename(hk_counts,c("Apiaceae.Umbell"="APIACEAE",
                                "Arceuthobium"="ARCEUTHOBI",
                                "Cyperaceae"="CYPERACE",
                                "Equisetum"="EQUISETU",
                                "Ericaceae"="ERICACEX",
                                "Euphorbia"="EUPHORB",
                                "Larix"="LARIXPSEU",
                                "Ostrya"="OSTRYCAR",
                                "Picea"="PICEAX",
                                "Polypodium"="POLYPOD",
                                "Ranunculus"="RANUNCUL",
                                "Rosaceae"="ROSACEAX",
                                "Rumex"="RUMEOXYR",
                                "Selaginella.rupestris"="SELAGINE",
                                "Tsuga"="TSUGAX",
                                "Urtica"="URTICACX"))  

ACERX = rowSums(hk_counts[,c("Acer.rubrum","Acer.saccharum","Acer.spicatum","Acer.undif.")]) #ACERX
ALNUSX = rowSums(hk_counts[,c("alnus.3.pore","Alnus.4","Alnus.5","AlnusUndif" )])#ALNUSX
BETULA = rowSums(hk_counts[,c("Betula","Betulaceae.undif")]) #"BETULA"     
FRAXINUX=rowSums(hk_counts[,c("Fraxinus3","Fraxinus4","FraxUnid")]) #"FRAXINUX"
IVA=rowSums(hk_counts[,c("Iva.ciliata","Iva.xanthifolia")]) #"IVA"   
JUGLANSX=rowSums(hk_counts[,c("Juglans.cinerea","Juglans.nigra","Juglans.undiff")])#"JUGLANSX"                             
LYCOPODX=rowSums(hk_counts[,c("Lyco.annotinum","Lyco.clavatum","Lyco.complanatum","Lyco.lucidulum","Lyco.obscurum","Lyco.selago","Lyco.undiff")])#"LYCOPODX" 
PINUSX=rowSums(hk_counts[,c("Pinus.subg..Pinus","Pinus.subg..Strobus","Pinus.undiff")])#"PINUSX"
POLYGONAX=rowSums(hk_counts[,c("Polygonum.amphibium.type","Polygonum.lapth..Type","Polygonum.sect..Persicaria")])#"POLYGONAX"    

cols_keep = which(toupper(colnames(hk_counts1))%in%colnames(pol.cal.count))
hk_counts2 = cbind(hk_counts[,1:7],hk_counts1[,cols_keep],ACERX,ALNUSX,
                   FRAXINUX,IVA,JUGLANSX,LYCOPODX,PINUSX,POLYGONAX)
hk_counts2[,"Betula"] = BETULA
colnames(hk_counts2) <- toupper(colnames(hk_counts2))
Other = rowSums(hk_counts[,8:ncol(hk_counts)]) - rowSums(hk_counts2[,8:ncol(hk_counts2)])

hk_counts3 = cbind(matrix(0,nrow(hk_counts2),10),hk_counts2,Other)

hk_counts3[hk_counts3$NAME=="DeepPine",12] <- c("Deep Pine")
hk_counts3$NAME<-sub(x=hk_counts3$NAME, pattern=c("DeerPrint"), replacement=c("Deer Print"))
hk_counts3$NAME<-sub(x=hk_counts3$NAME, pattern=c("CircleLily"), replacement=c("Circle Lily"))
hk_counts3$NAME<-sub(x=hk_counts3$NAME, pattern=c("Big Muskie07X"), replacement=c("BigMuskie07X"))
hk_counts3$NAME<-sub(x=hk_counts3$NAME, pattern=c("Trout07F"), replacement=c("Trout07X"))

#hk_counts3$NAME = as.factor(as.character(hk_counts3$NAME))

for(i in 1:nrow(hk_counts3)){
  hk_counts3[i,1:10] <- hk_meta[grep(hk_counts3$NAME[i],hk_meta$name),]
}

head(hk_counts3)

colnames(hk_counts3)<-c(colnames(hk_meta),colnames(hk_counts3[,11:ncol(hk_counts3)]))

hk_counts3 = as.data.frame(hk_counts3)
save(hk_counts3,file=paste0(data.dir,"hk_counts3.csv"))

load(paste0(data.dir,"hk_counts3.csv"))

hk_counts_bs<-compile_list_stepps(hk_counts3[,18:ncol(hk_counts3)], 'babystepps', pollen.equiv.stepps = pollen.equiv.stepps)

to_pol_mat = which(colnames(hk_counts_bs)%in%colnames(pol.cal.count))

add_on<-cbind(hk_counts3$name,hk_counts3$lat,hk_counts3$long,rep(0,length(hk_counts3$name)),hk_counts3$NAME,hk_counts3$AGE)
colnames(add_on) <- c("SiteID","LatitudeNorth","LongitudeWest","dataset.id","ContactName","Age")

pol.cal.count <- rbind(pol.cal.count,cbind(add_on,hk_counts_bs),deparse.level = 0)

##### Fix problems with matrix
rownames(pol.cal.count)<-seq(1,nrow(pol.cal.count),1)
pol.cal.count[is.na(pol.cal.count)]<-0
colnames(pol.cal.count)<-c("SiteID","LatitudeNorth","LongitudeWest","dataset.id","ContactName","Age",colnames(pol.cal.count[,7:ncol(pol.cal.count)]))
#save(pol.cal.count,file="pol.cal.count.mnwi1.csv")

#####
##### Plotting chronologies and maps of sites by age bins #####
#####

head(pol.cal.count)
ponds = pol.cal.count[pol.cal.count$Age>0,]
#ponds = ponds[ponds$Age<=5000,]
ponds = ponds[order(ponds$LatitudeNorth),]
site.factor = factor(ponds[,1],labels = seq(1,length(unique(ponds[,1])),1))
ponds1 = cbind(site.factor,ponds)

ponds_rm = which(ponds1$LatitudeNorth<44.5&ponds1$LongitudeWest>c(-86))
ponds1 = ponds1[-ponds_rm,]

pdf(paste0(fig.dir,"chrono.10k.pdf"))
par(mfrow=c(1,1))
#plot(ponds1$Age,ponds1$LatitudeNorth,pch=19,cex=.25,main="All Records",ylab="Latitude",xlab="Age BP")
#abline(v=2000)
plot(ponds1$Age,ponds1$LatitudeNorth,pch=19,cex=.5,xlim=c(0,10000),main="All Pollen Records",ylab="Latitude",xlab="Age BP")
dev.off()

plot.seq = seq(0,20000,500)

#quartz()
par(mfrow = c(2,2))
for(i in 2:length(plot.seq)){
  map('state', xlim=range(ponds1$LongitudeWest)+c(-2, 2), ylim=range(ponds1$LatitudeNorth)+c(-1, 1))
  points(ponds1[ponds1$Age>plot.seq[i-1]&ponds1$Age<plot.seq[i],]$LongitudeWest, 
         ponds1[ponds1$Age>plot.seq[i-1]&ponds1$Age<plot.seq[i],]$LatitudeNorth, 
         pch=19, cex=.5)
  title(c(plot.seq[i-1],"-",plot.seq[i]))
}
dev.off()




