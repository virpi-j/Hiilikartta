rm(list=ls())
gc()
if(length(dev.list())>0) dev.off()

toFile <- F

setwd("~/Hiilikartta")
nSitesRun0 <- 100
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

#source("~/HiilikarttaGH/settings.R")
devtools::source_url("https://raw.githubusercontent.com/virpi-j/Hiilikartta/master/settings.R")

#setX=1
#sampleIDs <- 1
set.seed(1)
#nn <- sample(1:dim(data.all)[1], nSitesRun, replace = F)
#dataS <- data.all[nn,]
nn0 <- order(data.all$age)[1:nSitesRun0] # The youngest segments to test set
nn1 <- sample(setdiff(1:dim(data.all)[1],nn0), nSitesRun - nSitesRun0)
dataS <- data.all[c(nn0,nn1),]
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
#source_url("https://raw.githubusercontent.com/virpi-j/Hiilikartta/master/functions.R")
source("~/HiilikarttaGH/functions.R")
outType <- "testRun"
harvScen <- "Base"
harvInten <- "Base"
manualRun <- F
if(manualRun){
  RCP=0
  harvScen <- "Base"
  harvInten <- "Base"
  easyInit=FALSE; forceSaveInitSoil=F; cons10run = F
  procDrPeat=F; coeffPeat1=-240; coeffPeat2=70
  coefCH4 = 0.34; coefN20_1 = 0.23; coefN20_2 = 0.077#g m-2 y-1
  landClassUnman=NULL; compHarvX = 0; funPreb = regionPrebas
  initSoilCreStart=NULL; outModReStart=NULL; reStartYear=1
  sampleX=NULL; P0currclim=NA; fT0=NA
  sampleID <- 1
}

print(paste("Run initialization for region",r_no,", sample size",nSitesRun))
out <- runModel(1,harvScen="Base",harvInten="Base",outType = outType)

CoeffSim <- T
ferti <- 1
fertmax <- 6
png("/scratch/project_2000994/PREBASruns/PREBAStesting/testPlots.png")
par(mfrow=c(ceiling(sqrt(fertmax)),floor(sqrt(fertmax))))
dataSorig <- dataS
output <- list()

for(ferti in 1:fertmax){
  print(paste("Fert =",ferti))
  toMem <- ls()
  dataS <- dataSorig
  if(CoeffSim){
    simInitData <- data.table(ba = 0.01, age = 1, dbh = 0.01, pine = 100, spruce = 0,
                              birch = 0, decid = 0, fert = ferti, h = 0.01, minpeat = 1,
                              landclass = 1, cons = 0)
    dataS$ba[1:nSitesRun0] <- simInitData$ba
    dataS$age[1:nSitesRun0] <- simInitData$age
    dataS$dbh[1:nSitesRun0] <- simInitData$dbh
    dataS$pine[1:nSitesRun0] <- simInitData$pine
    dataS$spruce[1:nSitesRun0] <- simInitData$spruce
    dataS$birch[1:nSitesRun0] <- simInitData$birch
    dataS$decid[1:nSitesRun0] <- simInitData$decid
    dataS$fert[1:nSitesRun0] <- simInitData$fert
    dataS$h[1:nSitesRun0] <- simInitData$h
    dataS$minpeat[1:nSitesRun0] <- simInitData$minpeat
    dataS$landclass[1:nSitesRun0] <- simInitData$landclass
    dataS$cons[1:nSitesRun0] <- simInitData$cons
    
    #simCol <- match(colnames(simInitData),colnames(dataS))
    #dataS[,c(simCol)] <- 
    #  data.table(matrix(simInitData,nSitesRun, ncol(simInitData), byrow = T))
    
    ops <- list(dataS)
  }
  #source_url("https://raw.githubusercontent.com/virpi-j/Hiilikartta/master/functions.R")
  source("~/HiilikarttaGH/functions.R")
  harvScen <- "baseTapio"
  harvScen = "Base"
  #harvScen <- "NoHarv"
  harvInten <- "Base"
  if(harvScen=="NoHarv") harvInten <- "NoHarv"
  if(manualRun){
    RCP=climScen; easyInit=FALSE; forceSaveInitSoil=F; cons10run = F; procDrPeat=F; coeffPeat1=-240; coeffPeat2=70; coefCH4 = 0.34; coefN20_1 = 0.23; coefN20_2 = 0.077; landClassUnman=NULL; compHarvX = 0; funPreb = regionPrebas; initSoilCreStart=NULL; outModReStart=NULL; reStartYear=1; sampleX=NULL; P0currclim=NA; fT0=NA; sampleID <- 1
  }
  
  if(!toFile){
    out <- lapply(sampleIDs, function(jx) {
      runModel(jx,harvScen=harvScen, harvInten=harvInten, outType = outType, RCP = climScen)})
  } else {
    out <- mclapply(sampleIDs, function(jx) {
      runModel(jx,harvScen=harvScen, harvInten=harvInten, outType = outType, RCP = climScen)
    }, mc.cores = nCores,mc.silent=FALSE)      
  }
  
  rm(list=setdiff(ls(),c(toMem,"out")))
  
  V <- age <- nep <- wTot <- wGV <- soilC <- array(0,c(nSitesRun0,nYears,length(sampleIDs)),
                                                   dimnames = list(paste0("site",1:nSitesRun0),
                                                                     2014+1:nYears,
                                                                     climMod))
  for(ij in 1:length(out)){
    #V[,,ij] <- apply(out[[ij]]$multiOut[,,"V",,1],1:2,"sum")
    V[,,ij] <- apply(out[[ij]]$region$multiOut[1:nSitesRun0,,"V",,1],1:2,"sum")
    age[,,ij] <- apply(out[[ij]]$region$multiOut[1:nSitesRun0,,"age",,1],1:2,"mean")
    nep[,,ij] <- apply(out[[ij]]$region$multiOut[1:nSitesRun0,,"NEP/SMI[layer_1]",,1],1:2,"sum")
    wTot[,,ij] <- apply(out[[ij]]$region$multiOut[1:nSitesRun0,,c(24,25,31,32,33),,1],1:2,"sum")
    wGV[,,ij] <- out[[ij]]$region$GVout[1:nSitesRun0,,4]
    soilC[,,ij] <- out[[ij]]$region$multiOut[1:nSitesRun0,,"soilC",1,1]
  }
  output[[ferti]] <- list(V, age, nep, wTot, wGV, soilC)
  names(output[[ferti]]) <- c("V", "age", "nep", "wTot", "wGV", "soilC")
  for(ij in 1:length(output[[ferti]])){
    if(names(output[[ferti]])[ij]!="age"){
      tmp <- output[[ferti]][[ij]]
      xi <- apply(tmp,2:3,mean)
      xmean <- apply(tmp,2,mean)
      agei <- apply(output[[ferti]]$age,2:3,mean)
      agemean <- apply(output[[ferti]]$age,2,mean)
      ylims <- c(0,max(xmean))
      plot(agemean,xmean,xlab="age",ylab=names(output[[ferti]])[ij],
           type="l",ylim=ylims,lwd=2,
           main=paste("Fert =",ferti))
      for(ik in sampleIDs){ # go through climate models
        #lines(as.numeric(names(Vmean)),Vi[,ij],col="blue")
        lines(agei[,ik],xi[,ik],col="blue")
      }
    }
  }
}
dev.off()
file.remove(paste0(path_initiSoilC,"initSoilCunc/forCent",r_no,"/initSoilC.rdata"))

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