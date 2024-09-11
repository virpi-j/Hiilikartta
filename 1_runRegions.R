rm(list=ls())
gc()
if(length(dev.list())>0) dev.off()

toFile <- F
set.seed(1)
setwd("~/Hiilikartta")
nSitesRun <-10000
nSitesRun0 <- 50
if(!toFile) nSitesRun <-1000
fertmax <- 6 # max fert type
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
startingYear = 2015
endingYear = 2050
nYears = endingYear-startingYear

#source_url("https://raw.githubusercontent.com/virpi-j/Hiilikartta/master/functions.R")
source("~/HiilikarttaGH/functions.R")
outType <- "testRun"
harvScen <- "Base"
harvInten <- "Base"
manualRun <- T
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
dataSorig <- dataS
output <- list()
startingYear = 2015
endingYear = 2100
nYears = endingYear-startingYear

species <- c(100,0,0,0) # pine, spruce, birch, deciduous
#harvScen <- "baseTapio"
#harvScen = "Base"
harvSceni <- "NoHarv"
harvScens <- c("NoHarv","Mitigation","BaseLow","adapt","baseTapio", "Base")
ferti <- 1
for(harvSceni in harvScens){
  for(initAge in c(NA,yearsToMem)){ # 50#NA
    for(ferti in 1:fertmax){
      time0 <- Sys.time()
      print(paste("fert =",ferti))
      toMem <- ls()
      dataS <- dataSorig
      if(CoeffSim){
        simInitData <- data.table(ba = 0.01, age = 1, dbh = 0.01, pine = species[1], 
                                  spruce =  species[2], birch =  species[3], decid =  species[4], 
                                  fert = ferti, h = 0.01, minpeat = 1,
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
      harvInten <- "Base"
      harvScen <- "Base"
      if(harvSceni!="Base"){ 
        harvScen <- harvSceni
        harvInten <- "Low"
      }
      if(harvSceni=="BaseLow"){ 
        harvScen <- "Base"
        harvInten <- "Low"
      }
      if(harvSceni=="NoHarv"){ 
        harvInten <- "NoHarv"
        harvScen <- "NoHarv"
      }
      if(manualRun){
        RCP=climScen; easyInit=FALSE; forceSaveInitSoil=F; cons10run = F; procDrPeat=F; coeffPeat1=-240; coeffPeat2=70; coefCH4 = 0.34; coefN20_1 = 0.23; coefN20_2 = 0.077; landClassUnman=NULL; compHarvX = 0; funPreb = regionPrebas; initSoilCreStart=NULL; outModReStart=NULL; reStartYear=1; sampleX=NULL; P0currclim=NA; fT0=NA; sampleID <- 1
      }
      if(!toFile){
        out <- lapply(sampleIDs, function(jx) {
          runModel(jx,harvScen=harvScen, harvInten=harvInten, outType = "hiiliKartta", RCP = climScen, initAge = initAge)})
      } else {
        out <- mclapply(sampleIDs, function(jx) {
          runModel(jx,harvScen=harvScen, harvInten=harvInten, outType = outType, RCP = climScen, initAge = initAge)
        }, mc.cores = nCores,mc.silent=FALSE)      
      }
      if(harvSceni=="NoHarv" & is.na(initAge)){
        multiOut <- array(0,dim = c(dim(out[[1]]$restartMod$multiOut),length(sampleIDs)))
        GVOut <- array(0,dim = c(dim(out[[1]]$restartMod$GVout),length(sampleIDs)))
        reStartSoil <- array(0,dim = c(dim(out[[1]]$reStartSoil),length(sampleIDs)))
        for(ij in sampleIDs){
          multiOut[,,,,,ij] <- out[[ij]]$restartMod$multiOut
          GVOut[,,,ij] <- out[[ij]]$restartMod$GVout
          reStartSoil[,,,,,ij] <- out[[ij]]$reStartSoil
        }
        multiOut <- apply(multiOut,1:(length(dim(multiOut))-1),mean)
        dimnames(multiOut) <- dimnames(out[[1]]$restartMod$multiOut)
        GVOut <- apply(GVOut,1:(length(dim(GVOut))-1),mean)
        dimnames(GVOut) <- dimnames(out[[1]]$restartMod$GVout)
        reStartSoil <- apply(reStartSoil,1:(length(dim(reStartSoil))-1),mean)
        dimnames(reStartSoil) <- dimnames(out[[1]]$reStartSoil)
        
        reStartMod <- list(multiOut, GVOut)
        names(reStartMod) <- c("multiOut", "GVOut")
        save(reStartMod,reStartSoil,
             file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/HiiliKartta_startStates",r_no,".rdata"))
        
      }
      rm(list=setdiff(ls(),c(toMem,"out")))
      
      V <- age <- nep <- wTot <- wGV <- soilC <- array(0,c(nSitesRun0,nYears,length(sampleIDs)),
                                                       dimnames = list(paste0("site",1:nSitesRun0),
                                                                       2014+1:nYears,
                                                                       climMod[sampleIDs]))
      for(ij in 1:length(out)){
        V[,,ij] <- out[[ij]]$V
        age[,,ij] <- out[[ij]]$age
        nep[,,ij] <- out[[ij]]$nep
        wTot[,,ij] <- out[[ij]]$wTot
        wGV[,,ij] <- out[[ij]]$wGV
        soilC[,,ij] <- out[[ij]]$soilC
      }
      output[[ferti]] <- list(V, age, nep, wTot, wGV, soilC)
      names(output[[ferti]]) <- c("V", "age", "nep", "wTot", "wGV", "soilC")
      print(Sys.time()-time0)
    }
    
    outFileePath <-"/scratch/project_2000994/PREBASruns/PREBAStesting/"
    if(is.na(initAge)) initAge <- 0
    outFilee <- paste0(outFileePath,"HiiliKarttaTestPlots_rno",r_no,"_",harvScen,"_",harvInten,
                       "_species",which(species>0),"_initAge",initAge,".pdf")
    pdf(outFilee)
    par(mfrow=c(ceiling(sqrt(fertmax)),floor(sqrt(fertmax))))
    for(ij in 1:length(output[[1]])){
      ymax <- max(output[[ferti]][[1]])
      for(ferti in 1:fertmax)   ymax <- max(ymax,max(colMeans(output[[ferti]][[ij]])))
      ymin <- min(output[[ferti]][[1]])
      for(ferti in 1:fertmax)   ymin <- min(ymin,min(colMeans(output[[ferti]][[ij]])))
      ylims <- c(ymin,ymax)
      
      for(ferti in 1:fertmax){
        if(names(output[[ferti]])[ij]!="age"){
          tmp <- output[[ferti]][[ij]]
          xi <- apply(tmp,2:3,mean)
          xmean <- apply(tmp,2,mean)
          timei <- 1:ncol(output[[ferti]][[ij]])
          agei <- apply(output[[ferti]]$age,2:3,mean)
          agemean <- apply(output[[ferti]]$age,2,mean)
          plot(agemean,xmean,xlab="age",ylab=names(output[[ferti]])[ij],
               type="l",ylim=ylims,lwd=2,
               main=paste("fert =",ferti))
          for(ik in sampleIDs){ # go through climate models
            #lines(as.numeric(names(Vmean)),Vi[,ij],col="blue")
            lines(agei[,ik],xi[,ik],col="blue")
          }
          if(ymin<0) lines(c(min(agei),max(agei)),c(0,0),col="black")
        }
      }
      
      for(ferti in 1:fertmax){
        tmp <- output[[ferti]][[ij]]
        xi <- apply(tmp,2:3,mean)
        xmean <- apply(tmp,2,mean)
        time <- 1:ncol(output[[ferti]][[ij]])
        #agei <- apply(output[[ferti]]$age,2:3,mean)
        #agemean <- apply(output[[ferti]]$age,2,mean)
        plot(time,xmean,xlab="time",ylab=names(output[[ferti]])[ij],
             type="l",ylim=ylims,lwd=2,
             main=paste("fert =",ferti))
        for(ik in sampleIDs){ # go through climate models
          #lines(as.numeric(names(Vmean)),Vi[,ij],col="blue")
          lines(time,xi[,ik],col="blue")
        }
        if(ymin<0) lines(c(min(time),max(time)),c(0,0),col="black")
      }
    }
    dev.off()
    outFilee <- paste0(outFileePath,"HiiliKarttaTestPlots_rno",r_no,"_",harvScen,"_",harvInten,
                       "_species",which(species>0),"_initAge",initAge,".rdata")
    save(output,file=outFilee)
  }
  
}  
###########

file.remove(paste0(path_initiSoilC,"initSoilCunc/forCent",r_no,"/initSoilC.rdata"))
file.remove(paste0(path_initiSoilC,"initSoilCunc/forCent",r_no,"/deadWV_mortMod",mortMod,".rdata"))

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