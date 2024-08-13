rm(list=ls())
gc()
if(length(dev.list())>0) dev.off()

toFile <- F

setwd("~/Hiilikartta")
if(!toFile) nSitesRun <-200
CSCrun <- T
vPREBAS <- "newVersion"
r_no <- 4
path_wrkdir  <- "/scratch/project_2000994/PREBASruns/finRuns/"
path_initiSoilC <- "/scratch/project_2000994/PREBASruns/finRuns/"
path_output <- "/scratch/project_2000994/PREBASruns/finRuns/"
landClassX <- 1:2

climatepath = "/scratch/project_2000994/RCP/"
climMod <- c("CanESM2.","CNRM.","GFDL.","HadGEM2.","MIROC.")
rcpx <- c("rcp26","rcp45","rcp85")

climScen <- 1 # 1=rcp2.6, 2=rcp4.5, 3=rcp8.5
climModid <- 1 # 
climModids <- sampleIDs <- 1:length(climMod) # for iterations
rcps <- rcpsFile <-paste0(climMod[climModid],rcpx[climScen])
rcpsName <- rcps

source("~/Hiilikartta/settings.R")
#devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/finRuns/Rsrc/settings.r")

#setX=1
#sampleIDs <- 1
set.seed(1)
nn <- sample(1:dim(data.all)[1], nSitesRun, replace = F)
#nSamples <- ceiling(dim(data.all)[1]/nSitesRun)
#sampleIDs <- split(1:nSamples,             # Applying split() function
#                   cut(seq_along(1:nSamples),
#                       nSetRuns,
#                       labels = FALSE))[[setX]]
dataS <- data.all[nn,]
#ops <- split(data.all, sample(1:nSamples, nrow(data.all), replace=T))
# test

load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/maakunta_",r_no,"_IDsTab.rdata"))
data.IDs <- data.IDs[segID!=0]
data.IDs$segID <- data.IDs$maakuntaID
setkey(data.IDs,segID)
setkey(dataS,segID)
#setkey(data.IDs,maakuntaID)

tabX <- merge(data.IDs,dataS)
ntabX <- tabX[,.I[which.max(y)],by=segID]$V1
#cbind(tabX[ntabX[1:10],"segID"],dataS[1:10,"segID"])
#tabX[,.I[which.max(tabX$y)],by=segID]$V1
#tabX[,.I[which.max(tabX$y)],by=segID]$V1
dataS <- cbind(dataS, tabX[ntabX,c("x","y")])

set_thin_PROJ6_warnings(TRUE)
xy <- dataS[,c("segID","x","y")]
coordinates(xy) <- c("x","y")
proj4string(xy) <- crsX
#cord = SpatialPoints(xy, proj4string=CRS("+init=EPSG:3067"))
location<-as.data.frame(spTransform(xy, CRS("+init=epsg:4326")))
dataS$lat <- location$y


ops <- list(dataS)
gc()
toMem <- ls()
###check and run missing sampleIDs 
# library('stringi')
# fileX <- list.files(path= "/scratch/project_2000994/PREBASruns/finRuns/outputDT/forCent12/", pattern = "age")
# sampleIDs <- which(!1:nSamples %in%  as.numeric(stri_extract_last(fileX, regex = "(\\d+)")))
# print(sampleIDs)
# sampleIDs <- c(66,342,395)
#source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/general/functions.r")
source("~/Hiilikartta/functions.R")
outType <- "testRun"
harvScen <- "Base"
harvInten <- "Base"

print(paste("Run initialization for region",r_no,", sample size",nSitesRun))
out <- runModel(1,harvScen="Base",harvInten="Base",outType = outType)

simInitData <- data.table(ba = 0.01, age = 1, dbh = 0.01, pine = 100, spruce = 0,
                          birch = 0, decid = 0, fert = 1, h = 0.01, minpeat = 1,
                          landclass = 1, cons = 0)

dataS[,match(colnames(simInitData),colnames(dataS))] <- 
  data.table(matrix(simInitData,nSitesRun, ncol(simInitData), byrow = T))

ops <- list(dataS)
source("~/Hiilikartta/functions.R")
harvScen <- "baseTapio"
harvInten <- "Base"

if(!toFile){
  out <- lapply(sampleIDs, function(jx) {
    runModel(jx,harvScen=harvScen, harvInten=harvInten, outType = outType, RCP = climScen)})
} else {
  out <- mclapply(sampleIDs, function(jx) {
    runModel(jx,harvScen=harvScen, harvInten=harvInten, outType = outType, RCP = climScen)
  }, mc.cores = nCores,mc.silent=FALSE)      
}

rm(list=setdiff(ls(),c(toMem,"out")))
file.remove(paste0(path_initiSoilC,"initSoilCunc/forCent",r_no,"/initSoilC.rdata"))

V <- array(0,c(nSitesRun,nYears,length(sampleIDs)),dimnames = list(paste0("site",1:nSitesRun),
                                                                   2014+1:nYears,
                                                                   climMod))
for(ij in 1:length(out)){
  V[,,ij] <- apply(out[[ij]]$region$multiOut[,,"V",,1],1:2,"sum")
}
Vi <- apply(V,2:3,mean)
Vmean <- apply(V,2,mean)
plot(as.numeric(names(Vmean)),Vmean,type="l",ylim=c(min(Vi),max(Vi)),lwd=2)
for(ij in sampleIDs){
  lines(as.numeric(names(Vmean)),Vi[,ij],col="blue")
}

# models outputs to NAs, outputDT, initSoilC and plots
Sys.chmod(list.dirs("NAs"), "0777",use_umask=FALSE)
f <- list.files("NAs", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)

Sys.chmod(list.dirs("outputDT"), "0777",use_umask=FALSE)
f <- list.files("outputDT", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)

Sys.chmod(list.dirs("initSoilC"), "0777",use_umask=FALSE)
f <- list.files("initSoilC", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)

Sys.chmod(list.dirs("plots"), "0777",use_umask=FALSE)
f <- list.files("plots", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)