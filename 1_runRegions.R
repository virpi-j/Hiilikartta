#rm(list=ls())
#gc()
if(length(dev.list())>0) dev.off()
print(paste("region",r_no))
toFile <- F
set.seed(1)
setwd("~/HiilikarttaGH/")
#source("/scratch/project_2000994/PREBASruns/PREBAStesting/localSettins.R",local=T)

nSitesRun <-100#00
nSitesRun0 <- 10#0
fertmax <- 5 # max fert type
if(testaus){
  nSitesRun <-100
  nSitesRun0 <- 10
  fertmax <- 2 # max fert type
  yearsToMem <- c(30,50)
}
CSCrun <- T
#vPREBAS <- "newVersion"
#r_no <- 4
path_wrkdir  <- "/scratch/project_2000994/PREBASruns/finRuns/"
path_initiSoilC <- "/scratch/project_2000994/PREBASruns/finRuns/"
path_output <- "/scratch/project_2000994/PREBASruns/finRuns/"
landClassX <- 1:2

climatepath = "/scratch/project_2000994/RCP/"
climMod <- c("CanESM2.","CNRM.","GFDL.","HadGEM2.","MIROC.")
rcpx <- c("rcp26","rcp45","rcp85")
if(testaus){
  climMod <- c("CanESM2.","CNRM.")
}
climScen <- 1 # 1=rcp2.6, 2=rcp4.5, 3=rcp8.5
climModid <- 1 # 
climModids <- sampleIDs <- 1:length(climMod) # for iterations
rcps <- rcpsFile <-paste0(climMod[climModid],rcpx[climScen])
rcpsName <- rcps

source("~/HiilikarttaGH/settings.R", local = T)
print(paste("vPREBAS =", vPREBAS))
regnames <- c("Uusimaa", "Ahvenanmaa", "Keski-Pohjanmaa", "Pirkanmaa",
              "Etela-Karjala", "Keski-Suomi", "Pohjois-Savo", 
              "Lappi", "Kanta-Hame", "Pohjanmaa", "Varsinais-Suomi",
              "Etela-pohjanmaa", "Paijat-Hame", "Satakunta", "Kymenlaakso",
              "Kainuu", "Etela-Savo", "Pohjois-Karjala", "Pohjois-Pohjanmaa")
r_nos <- c(1, 21, 16, 6, 9, 13, 11, 19, 5, 15, 2, 14, 7, 4, 8, 18, 10, 12, 17) # Official numbering compared to PREBAS run region ids
#devtools::source_url("https://raw.githubusercontent.com/virpi-j/Hiilikartta/master/settings.R")

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
endingYear = 2100
nYears = endingYear-startingYear

#source_url("https://raw.githubusercontent.com/virpi-j/Hiilikartta/master/functions.R")
source("~/HiilikarttaGH/functions.R", local = T)
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
  sampleID = 1; initAge=NA
}

# initialize 
print(paste("Run initialization for region",r_no,", sample size",nSitesRun))
out <- runModel(1,harvScen="Base",harvInten="Base",outType = outType)

CoeffSim <- T
ferti <- 1
dataSorig <- dataS
output <- list()
startingYear = 2015
endingYear = 2100
nYears = endingYear-startingYear

ages <- data.all[,"age"]
tmps <- data.all[,c("pine","spruce","birch","fert","age","ba")]
tmp <- tmps[,1:3]/rowSums(tmps[,1:3])
tmp[rowSums(tmps)==0,]<-0
nsets <- 7
speciess <- array(0,c(nSitesRun0,3,nsets,fertmax),
                dimnames = list(1:nSitesRun0,c("pine","spruce","birch"),  
                                c("pinedom","sprucedom","birchdom",
                                                "pinebirch","sprucebirch","sprucepine","sprucepinebirch"),
                                paste0("fert",1:fertmax)))
for(ferti in 1:fertmax){
  # pine dominated
  pinedom <- tmps[which(tmp$pine>0.5 & tmps$fert==ferti),]
  pinedom <- pinedom[sample(1:nrow(pinedom),size = nSitesRun0, replace = T),1:3]
  speciess[,,"pinedom",ferti] <- as.matrix(pinedom)
  # spruce dominated
  sprucedom <- tmps[which(tmp$spruce>0.5 & tmps$fert==ferti),]
  sprucedom <- sprucedom[sample(1:nrow(sprucedom),size = nSitesRun0, replace = T),1:3]
  speciess[,,"sprucedom",ferti] <- as.matrix(sprucedom)
  # birch dominated
  birchdom <- tmps[which(tmp$birch>0.5 & tmps$fert==ferti),]
  birchdom <- birchdom[sample(1:nrow(birchdom),size = nSitesRun0, replace = T),1:3]
  speciess[,,"birchdom",ferti] <- as.matrix(birchdom)
  # pine & birch mix
  pinebirch <- tmps[which(tmp$pine<0.5 & tmp$birch<0.5 & tmp$spruce<0.2 & tmps$fert==ferti),]
  pinebirch <- pinebirch[sample(1:nrow(pinebirch),size = nSitesRun0, replace = T),1:3]
  speciess[,,"pinebirch",ferti] <- as.matrix(pinebirch)
  # spruce & birch mix
  sprucebirch <- tmps[which(tmp$spruce<0.5 & tmp$birch<0.5 & tmp$pine<0.2 & tmps$fert==ferti),]
  sprucebirch <- sprucebirch[sample(1:nrow(sprucebirch),size = nSitesRun0, replace = T),1:3]
  speciess[,,"sprucebirch",ferti] <- as.matrix(sprucebirch)
  # spruce & pine mix
  sprucepine <- tmps[which(tmp$spruce<0.5 & tmp$pine<0.5 & tmp$birch<0.2 & tmps$fert==ferti),]
  sprucepine <- sprucepine[sample(1:nrow(sprucepine),size = nSitesRun0, replace = T),1:3]
  speciess[,,"sprucepine",ferti] <- as.matrix(sprucepine)
  # all mix
  sprucepinebirch <- tmps[which(tmp$spruce<.5 & tmp$pine<.5 & tmp$birch<.5 & tmp$spruce>0.2 & tmp$pine>.2 & tmp$birch>0.2 & tmps$fert==ferti),]
  sprucepinebirch <- sprucepinebirch[sample(1:nrow(sprucepinebirch),size = nSitesRun0, replace = T),1:3]
  speciess[,,"sprucepinebirch",ferti] <- as.matrix(sprucepinebirch)
}

#speciess <- array(0,c(nsets,4),dimnames = list(paste0("iter",1:nsets),c("pine","spruce","birch","decid")))
#speciess[1,] <- c(100,0,0,0)#tmps[which.max(tmp$pine),]#c(100,5,5,0)
#speciess[2,] <- c(0,100,0,0)#tmps[which.max(tmp$spruce),]#c(5,100,5,0) # pine, spruce, birch, deciduous
#speciess[3,] <- c(0,0,100,0)#tmps[which.max(tmp$birch),]#c(5,5,100,0)
#speciess[4,] <- c(0,50,150,0)#c(2,49,49,0)
#speciess[5,] <- c(50,50,0,0)#c(49,2,49,0)
speciesNamesLong <- speciesNames <- dimnames(speciess)[[3]]#c("100sp1","100sp2","sp3","50sp250sp3","50sp150sp3")
#speciesNamesLong <- c("manty", "kuusi", "lehtipuu", "kuusi-lehtipuu","manty-lehtipuu")
landclass0 <- 1 # landclass for sample0
minpeat0 <- 1 # mineral soils for sample0 (minral, drainet peat, undrained peat)
soiltype0 <- 1 # soiltypes (mineral soils, spruce mire, pine mire, ombrotrophic bog)

harvSceni <- "NoHarv"
harvScens <- c("NoHarv","Mitigation","BaseLow","adapt","baseTapio", "Base")
harvScens <- c("NoHarv","baseTapio")
#harvScens <- c("baseTapio","NoHarv")
ferti <- 1
speciesSeti <- 1
speciesName <- speciesNames[speciesSeti]
runPerHarvScen <- function(harvSceni, speciesSeti, dataS=dataSorig){
  if(harvSceni=="Base"){ 
    harvInten <- "Base"
    harvScen <- "Base"
  } else if(harvSceni=="BaseLow"){ 
    harvScen <- "Base"
    harvInten <- "Low"
  } else if(harvSceni=="NoHarv"){ 
    harvInten <- "NoHarv"
    harvScen <- "NoHarv"
  } else { 
    harvScen <- harvSceni
    harvInten <- "Low"
  }
  n0only <- F
  if(harvSceni%in%c("NoHarv","baseTapio")) n0only <- T
  print(paste("Species",speciesName,"run..."))
  for(ferti in 1:fertmax){
    ferti <<- ferti # make global
    print(paste(harvScen,"/",harvInten,"/ fert =",ferti))
    time0 <- Sys.time()
    outputAgei <-list()
    init0 <- 1
    if(ingrowth) init0 <- 0
    initAgei <- init0
    for(initAgei in init0:length(c(NA,yearsToMem))){ # 50#NA
      if(initAgei==0){ 
        initAge <- NA
      } else {initAge <- c(NA,yearsToMem)[initAgei]}
      toMem <- ls()
      dataS <- dataSorig
      if(CoeffSim){
        simInitData <- data.table(ba = 0.01, age = 0, dbh = 0.01, pine = 1, 
                                  spruce =  1, birch =  1, #decid =  species[4], 
                                  fert = ferti, h = 1.35, minpeat = minpeat0,
                                  landclass = landclass0, cons = 0)
        dataS$ba[1:nSitesRun0] <- simInitData$ba
        dataS$age[1:nSitesRun0] <- simInitData$age
        dataS$dbh[1:nSitesRun0] <- simInitData$dbh
        dataS$pine[1:nSitesRun0] <- speciess[,"pine",speciesSeti,ferti]#simInitData$pine
        dataS$spruce[1:nSitesRun0] <- speciess[,"spruce",speciesSeti,ferti]#simInitData$spruce
        dataS$birch[1:nSitesRun0] <-speciess[,"birch",speciesSeti,ferti]# simInitData$birch
        dataS$decid[1:nSitesRun0] <- speciess[,"birch",speciesSeti,ferti]#simInitData$decid
        dataS$fert[1:nSitesRun0] <- simInitData$fert
        dataS$h[1:nSitesRun0] <- simInitData$h
        dataS$minpeat[1:nSitesRun0] <- simInitData$minpeat
        dataS$landclass[1:nSitesRun0] <- simInitData$landclass
        dataS$cons[1:nSitesRun0] <- simInitData$cons
        if(n0only) dataS <- dataS[1:nSitesRun0,]
        ops <<- list(dataS)
      }
      
      #source_url("https://raw.githubusercontent.com/virpi-j/Hiilikartta/master/functions.R")
      #source("~/HiilikarttaGH/functions.R", local = T)
      if(manualRun){
        RCP=climScen; easyInit=FALSE; forceSaveInitSoil=F; cons10run = F; procDrPeat=F; outType = "hiiliKartta"; coeffPeat1=-240; coeffPeat2=70; coefCH4 = 0.34; coefN20_1 = 0.23; coefN20_2 = 0.077; landClassUnman=NULL; compHarvX = 0; funPreb = regionPrebas; initSoilCreStart=NULL; outModReStart=NULL; reStartYear=1; sampleX=NULL; P0currclim=NA; fT0=NA; sampleID <- 1
      }
      if(initAgei==0){
        out <- lapply(sampleIDs, function(jx) {
          runModel(jx,harvScen=harvScen, harvInten=harvInten, outType = "hiiliKartta", 
                   RCP = climScen, initAge = initAge)})
        #} else {
        #  out <- mclapply(sampleIDs, function(jx) {
        #    runModel(jx,harvScen=harvScen, harvInten=harvInten, outType = outType, RCP = climScen, initAge = initAge)
        #  }, mc.cores = nCores,mc.silent=FALSE)      
        #}
        if(harvSceni==harvScens[1] & initAgei==0){ # is.na(initAge)){
          multiOut <- array(0,dim = c(dim(out[[1]]$restartMod$multiOut),length(sampleIDs)))
          GVOut <- array(0,dim = c(dim(out[[1]]$restartMod$GVout),length(sampleIDs)))
          reStartSoil <- array(0,dim = c(dim(out[[1]]$reStartSoil),length(sampleIDs)))
          for(ij in sampleIDs){
            if(length(dim(multiOut))==5){
              multiOut[,,,,ij] <- out[[ij]]$restartMod$multiOut
              GVOut[,,,ij] <- out[[ij]]$restartMod$GVout
              reStartSoil[,,,,ij] <- out[[ij]]$reStartSoil
            } else {
              multiOut[,,,,,ij] <- out[[ij]]$restartMod$multiOut
              GVOut[,,,ij] <- out[[ij]]$restartMod$GVout
              reStartSoil[,,,,,ij] <- out[[ij]]$reStartSoil
            }
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
               file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/HiiliKartta_startStates",r_no,"_fert",ferti,".rdata"))
          
        }
      #}
      } else {
        out <- lapply(sampleIDs, function(jx) {
          runModel(jx,harvScen=harvScen, harvInten=harvInten, outType = "hiiliKartta", ingrowth = ingrowth, 
                   RCP = climScen, initAge = initAge)})
      }
      
      rm(list=setdiff(ls(),c(toMem,"out")))
      
      if(initAgei>0){
        V <- litters <- age <- nep <- wTot <- wGV <- soilC <- Vpine <- Vspruce <- Vbirch <-
          array(0,c(nSitesRun0,nYears,length(sampleIDs)),
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
          litters[,,ij] <- out[[ij]]$litters
          Vpine[,,ij] <- out[[ij]]$Vpine
          Vspruce[,,ij] <- out[[ij]]$Vspruce
          Vbirch[,,ij] <- out[[ij]]$Vbirch
        }
        outputAgei[[initAgei]] <- list(V, age, nep, wTot, wGV, soilC, litters, Vpine, Vspruce, Vbirch)
        names(outputAgei[[initAgei]]) <- c("V", "age", "nep", "wTot", "wGV", "soilC", "litters", "Vpine", "Vspruce", "Vbirch")
        print(Sys.time()-time0)
        plot(V[1,,1],type="l", ylim = c(-0.5, max(V)))
        lines(Vpine[1,,1],col="red")
        lines(Vspruce[1,,1],col="blue")
        lines(Vbirch[1,,1],col="green")
      }
    }
    output[[ferti]] <- outputAgei
    names(output[[ferti]]) <- paste0("initAge",c(0, yearsToMem))
    print(paste0(harvScen," / fert",ferti," / age0 / V:"))
    print(output[[ferti]][[1]]$V[1,1:10,1])
    print(paste0(harvScen," / fert",ferti," / age30 / V:"))
    print(output[[ferti]][[2]]$V[1,1:10,1])
  }
  names(output) <- paste0("fert",1:fertmax)
  
  outFileePath <-"/scratch/project_2000994/PREBASruns/PREBAStesting/HiilikarttaResults/"
  #if(is.na(initAge)) initAge <- 0
  outFilee <- paste0(outFileePath,"HiiliKarttaTestPlots_rno",r_no,"_",harvScen,"_",harvInten,
                     "_species",speciesName,".pdf")
  pdf(outFilee)
  for(ij in 1:length(output[[1]][[1]])){
    par(mfrow=c(ceiling(sqrt(fertmax)),floor(sqrt(fertmax))))
    ymax <- max(output[[1]][[1]][[ij]])
    ymin <- min(output[[1]][[1]][[ij]])
    for(ferti in 1:fertmax){
      for(agei in 1:length(c(0,yearsToMem))){
        ymax <- max(ymax,max(colMeans(output[[ferti]][[agei]][[ij]])))
        ymin <- min(ymin,min(colMeans(output[[ferti]][[agei]][[ij]])))
      }
    }
    ylims <- c(ymin,ymax)
   
    ageplot <- F
    if(ageplot){
    for(agei in 1:length(c(0,yearsToMem))){
      for(ferti in 1:fertmax){
        if(names(output[[ferti]][[agei]])[ij]!="age"){
          tmp <- output[[ferti]][[agei]][[ij]]
          xi <- apply(tmp,2:3,mean)
          xmean <- apply(tmp,2,mean)
          timei <- 1:ncol(output[[ferti]][[agei]][[ij]])
          ages <- apply(output[[ferti]][[agei]]$age,2:3,mean)
          agemean <- apply(output[[ferti]][[agei]]$age,2,mean)
          plot(agemean,xmean,xlab="age",ylab=names(output[[ferti]][[agei]])[ij],
               type="l",ylim=ylims,lwd=2,
               main=paste("fert =",ferti,"init year",c(0,yearsToMem)[agei]))
          for(ik in sampleIDs){ # go through climate models
            #lines(as.numeric(names(Vmean)),Vi[,ij],col="blue")
            lines(ages[,ik],xi[,ik],col="blue")
          }
          if(ymin<0) lines(c(min(ages),max(ages)),c(0,0),col="black")
        }
      }
    }
      }
    for(agei in 1:length(c(0,yearsToMem))){
      for(ferti in 1:fertmax){
        tmp <- output[[ferti]][[agei]][[ij]]
        xi <- apply(tmp,2:3,mean)
        xmean <- apply(tmp,2,mean)
        time <- 1:ncol(output[[ferti]][[agei]][[ij]])
        plot(time,xmean,xlab="time",ylab=names(output[[ferti]][[agei]])[ij],
             type="l",ylim=ylims,lwd=2,
             main=paste("fert =",ferti,"init year",c(0,yearsToMem)[agei]))
        for(ik in sampleIDs){ # go through climate models
          #lines(as.numeric(names(Vmean)),Vi[,ij],col="blue")
          lines(time,xi[,ik],col="blue")
        }
        if(ymin<0) lines(c(min(time),max(time)),c(0,0),col="black")
      }
    }
  }
  dev.off()
  outFilee <- paste0(outFileePath,"HiiliKarttaTestPlots_rno",r_no,"_",harvScen,"_",harvInten,
                     "_species",speciesName,".rdata")
  save(output,file=outFilee)
  print(paste("Saved file",outFilee))
}
###########
ij <- 1

for(ij in 1:nrow(speciess)){
  #species <<- speciess[ij,]
  speciesName <<- speciesNames[ij]
  runOut <- lapply(harvScens[1], function(jx) {
    runPerHarvScen(jx,speciesSeti = ij)})
  if(!parRuns){
    runOut <- lapply(harvScens[-1], function(jx) {
      runPerHarvScen(jx,speciesSeti = ij)})
  } else {
    runOut <- mclapply(harvScens[-1], function(jx) {
      runPerHarvScen(jx,speciesSeti = ij)
    }, mc.cores = 5, mc.silent=FALSE)      
  }      
}
  


###########
if(FALSE){
  ijs <- 1
  dfs <- data.frame()
  for(harvSceni in harvScens){
    if(harvSceni=="Base"){ 
      harvInten <- "Base"
      harvScen <- "Base"
    } else if(harvSceni=="BaseLow"){ 
      harvScen <- "Base"
      harvInten <- "Low"
    } else if(harvSceni=="NoHarv"){ 
      harvInten <- "NoHarv"
      harvScen <- "NoHarv"
    } else { 
      harvScen <- harvSceni
      harvInten <- "Low"
    }
    for(ijs in 1:nrow(speciess)){
      species <<- speciess[ijs,]
      speciesName <<- speciesNames[ijs]
      
      outFileePath <-"/scratch/project_2000994/PREBASruns/PREBAStesting/"
      outFilee <- paste0(outFileePath,"HiiliKarttaTestPlots_rno",r_no,"_",harvScen,"_",harvInten,
                         "_species",speciesName,".rdata")
      load(file=outFilee)
      
      # datatables...
      for(ferti in 1:fertmax){
        age <- round(apply(output[[1]]$initAge0$age[,,],2,mean))
        tmp <- apply(output[[ferti]]$initAge0$wTot,2,mean)/1000
        df <- data.table(Region=r_nos[r_no], Regname = regnames[r_no],
                         Maingroup = landclass0, Soiltype = soiltype0,
                         Drainage = minpeat0, Fertility = ferti, 
                         Species = ijs, SpeciesNames = speciesNamesLong[ijs],
                         Structure = 1, Regime = harvSceni, 
                         Rotation = 0)
        df2 <- data.table(t(tmp))
        colnames(df2) <- as.character(1:length(tmp))#as.character(age)
        df <- cbind(df, df2)
        dfs <- rbind(dfs, df)
      }
      
    
      #if(is.na(initAge)) initAge <- 0
      outFilee <- paste0(outFileePath,"HiiliKarttaTestPlots_rno",r_no,"_",harvScen,"_",harvInten,
                       "_species",which(species>0),".pdf")
      pdf(outFilee)
      par(mfrow=c(ceiling(sqrt(fertmax)),floor(sqrt(fertmax))))
      for(ij in 1:length(output[[1]][[1]])){
        ymax <- max(output[[1]][[1]][[ij]])
        ymin <- min(output[[1]][[1]][[ij]])
        for(ferti in 1:fertmax){
          for(agei in 1:length(c(0,yearsToMem))){
            ymax <- max(ymax,max(colMeans(output[[ferti]][[agei]][[ij]])))
            ymin <- min(ymin,min(colMeans(output[[ferti]][[agei]][[ij]])))
          }
        }
        ylims <- c(ymin,ymax)
        
        for(agei in 1:length(c(0,yearsToMem))){
          for(ferti in 1:fertmax){
            if(names(output[[ferti]][[agei]])[ij]!="age"){
              tmp <- output[[ferti]][[agei]][[ij]]
              xi <- apply(tmp,2:3,mean)
              xmean <- apply(tmp,2,mean)
              timei <- 1:ncol(output[[ferti]][[agei]][[ij]])
              ages <- apply(output[[ferti]][[agei]]$age,2:3,mean)
              agemean <- apply(output[[ferti]][[agei]]$age,2,mean)
              plot(agemean,xmean,xlab="age",ylab=names(output[[ferti]][[agei]])[ij],
                   type="l",ylim=ylims,lwd=2,
                   main=paste("fert =",ferti,"init year",c(0,yearsToMem)[agei]))
              for(ik in sampleIDs){ # go through climate models
                #lines(as.numeric(names(Vmean)),Vi[,ij],col="blue")
                lines(ages[,ik],xi[,ik],col="blue")
              }
              if(ymin<0) lines(c(min(ages),max(ages)),c(0,0),col="black")
            }
          }
        }
        for(agei in 1:length(c(0,yearsToMem))){
          for(ferti in 1:fertmax){
            tmp <- output[[ferti]][[agei]][[ij]]
            xi <- apply(tmp,2:3,mean)
            xmean <- apply(tmp,2,mean)
            time <- 1:ncol(output[[ferti]][[agei]][[ij]])
            plot(time,xmean,xlab="time",ylab=names(output[[ferti]][[agei]])[ij],
                 type="l",ylim=ylims,lwd=2,
                 main=paste("fert =",ferti,"init year",c(0,yearsToMem)[agei]))
            for(ik in sampleIDs){ # go through climate models
              #lines(as.numeric(names(Vmean)),Vi[,ij],col="blue")
              lines(time,xi[,ik],col="blue")
            }
            if(ymin<0) lines(c(min(time),max(time)),c(0,0),col="black")
          }
        }
      }
      dev.off()
    }    
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
