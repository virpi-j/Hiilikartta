#rm(list=ls())
#gc()
if(length(dev.list())>0) dev.off()
print(paste("region",r_no))
set.seed(1)
setwd(projDir)
#source("/scratch/project_2000994/PREBASruns/PREBAStesting/localSettins.R",local=T)

if(!exists("nSitesRun")) nSitesRun <-10000
nSitesRun0 <- 100
fertmax <- fertmax0 <- 5 # max fert type
if(testaus){
  nSitesRun <-2000
  nSitesRun0 <- 100
  fertmax <- fertmax0 <- 2 # 5 max fert type
  yearsToMem <- c(30,50) #c(30,50,70)
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

source(paste0(projDir,"settings.R"))#, local = T)
print(paste("vPREBAS =", vPREBAS))
regnames <- c("Uusimaa", "Ahvenanmaa", "Keski-Pohjanmaa", "Pirkanmaa",
              "Etela-Karjala", "Keski-Suomi", "Pohjois-Savo", 
              "Lappi", "Kanta-Hame", "Pohjanmaa", "Varsinais-Suomi",
              "Etela-pohjanmaa", "Paijat-Hame", "Satakunta", "Kymenlaakso",
              "Kainuu", "Etela-Savo", "Pohjois-Karjala", "Pohjois-Pohjanmaa")
r_nos <- c(1, 21, 16, 6, 9, 13, 11, 19, 5, 15, 2, 14, 7, 4, 8, 18, 10, 12, 17) # Official numbering compared to PREBAS run region ids
#devtools::source_url("https://raw.githubusercontent.com/virpi-j/Hiilikartta/master/settings.R")

set.seed(1)
nn0 <- order(data.all$age)[1:nSitesRun0] # The youngest segments to test set
nn1 <- sample(setdiff(1:dim(data.all)[1],nn0), nSitesRun - nSitesRun0)
dataS <- data.all[c(nn0,nn1),]

print(paste("Region",regnames[r_no]))
load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/maakunta_",r_no,"_IDsTab.rdata"))
data.IDs <- data.IDs[segID!=0]
data.IDs$segID <- data.IDs$maakuntaID
setkey(data.IDs,segID)
setkey(dataS,segID)

tabX <- merge(data.IDs,dataS)
ntabX <- tabX[,.I[which.max(y)],by=segID]$V1
dataS <- cbind(dataS, tabX[ntabX,c("x","y")])

#set_thin_PROJ6_warnings(TRUE)
xy <- dataS[,c("segID","x","y")]
coordinates(xy) <- c("x","y")
proj4string(xy) <- crsX
location<-as.data.frame(spTransform(xy, CRS("+init=epsg:4326")))
dataS$lat <- location$coords.x2

#ops <- list(dataS)
gc()
toMem <- ls()
startingYear = 2015
endingYear = 2100
nYears = endingYear-startingYear

#source_url("https://raw.githubusercontent.com/virpi-j/Hiilikartta/master/functions.R")
source(paste0(projDir,"functions.R"), local = T)
outType <- "testRun"
harvScen <- "Base"
harvInten <- "Base"
manualRun <- F
if(manualRun){
  RCP=0
  harvScen <- "Base"
  harvInten <- "Base"
  rcps<-"CurrClim"
  easyInit=FALSE; forceSaveInitSoil=F; cons10run = F
  procDrPeat=F; coeffPeat1=-240; coeffPeat2=70
  coefCH4 = 0.34; coefN20_1 = 0.23; coefN20_2 = 0.077#g m-2 y-1
  landClassUnman=NULL; compHarvX = 0; funPreb = regionPrebas
  initSoilCreStart=NULL; outModReStart=NULL; reStartYear=1
  sampleX = dataS; P0currclim=NA; fT0=NA
  climdata=NULL
  sampleID = 1; initAge=NA
  sampleX <- dataS
}

# initialize 
print(paste("Run initialization for region",r_no,", sample size",nSitesRun,"..."))
out <- runModel(1,sampleX = dataS, harvScen="Base",harvInten="Base",rcps="CurrClim",RCP=0, outType = outType)
print("done.")
clim <- out$clim
#colMeans(apply(out$region$multiOut[which(dataS$landclass==2),1:10,"grossGrowth",,1],1:2,sum))
#colMeans(apply(out$region$multiOut[which(dataS$landclass==1),1:10,"soilC",,1],1:2,sum))
#colMeans(apply(out$region$multiOut[which(dataS$landclass==2),1:10,"soilC",,1],1:2,sum))
CoeffSim <- T
ferti <- 1
dataSorig <- dataS
output <- list()
startingYear = 2015
endingYear = 2100
nYears = endingYear-startingYear

ages <- data.all[,"age"]
tmps <- data.all[,c("pine","spruce","birch","fert","age","ba","landclass")]
tmp <- tmps[,1:3]/rowSums(tmps[,1:3])
tmp[rowSums(tmps)==0,]<-0
#tmps[,1:3] <- tmp
nsets <- 9

speciess <- array(0,c(nSitesRun0,3,nsets,fertmax),
                  dimnames = list(1:nSitesRun0,c("pine","spruce","birch"),  
                                  c("typical","kitu","pinedom","sprucedom","birchdom",
                                    "pinebirch","sprucebirch","sprucepine","sprucepinebirch"),
                                  paste0("fert",1:fertmax)))

# typical forest in region
typical <- tmps[sample(1:nrow(tmps),size = nSitesRun0, replace = T),]
ferttypical <- typical$fert

# typical kitumaa/poorly productive forest in region
kitu <- tmps[sample(which(data.all$landclass==2),nSitesRun0, replace = T),]
kitu$fert[which(kitu$fert<5)] <- 5
fertkitu <- kitu$fert

for(ferti in 1:fertmax){
  # typical forest in region
  speciess[,,"typical",ferti] <- as.matrix(typical[,1:3])
  # typical forest in region
  speciess[,,"kitu",ferti] <- as.matrix(kitu[,1:3])
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

speciesNamesLong <- speciesNames <- dimnames(speciess)[[3]]#c("100sp1","100sp2","sp3","50sp250sp3","50sp150sp3")
landclass0 <- 1 # landclass for sample0
minpeat0 <- 1 # mineral soils for sample0 (minral, drainet peat, undrained peat)
soiltype0 <- 1 # soiltypes (mineral soils, spruce mire, pine mire, ombrotrophic bog)

harvSceni <- harvScens[3] #"NoHarv","Recreation","baseTapio","powerline_under","powerline_border"
ferti <- 1
speciesSeti <- 6
speciesName <- speciesNames[speciesSeti]
runPerHarvScen <- function(harvSceni, speciesSeti, dataS=dataSorig){
  if((speciesSeti==2 & harvSceni=="NoHarv") | speciesSeti!=2){
    ingrowth <- T
    if(harvSceni=="Base"){ 
      harvInten <- "Base"
      harvScen <- "Base"
    } else if(harvSceni=="BaseLow"){ 
      harvScen <- "Base"
      harvInten <- "Low"
    } else if(harvSceni%in%c("NoHarv")){ 
      harvInten <- "NoHarv"
      harvScen <- "NoHarv"
      #if(speciesName=="kitu"){
      #  initAges <- NA
      #}
    } else if (harvSceni%in%c("Powerline_under","Powerline_border")) { 
      harvScen <- harvSceni
      harvInten <- "Low"
      initAges <- NA
      ingrowth <- F
    } else { 
      harvScen <- harvSceni
      harvInten <- "Low"
    }
    print(paste(harvScen,"/",harvInten))
    n0only <- F # T if run only modified segments (no dependence on harvest rate in other segments)
    if(harvSceni%in%c("NoHarv","Recreation","baseTapio","Powerline_under","Powerline_border")) n0only <- T
    print(paste("Species",speciesName,"run..."))
    inAs <- c(0,yearsToMem)
    if(harvSceni%in%c("Powerline_under","Powerline_border")) inAs <- 0
    #if(speciesName=="kitu") inAs <- 0
    fertmax <- fertmax0
    if(speciesName%in%c("typical","kitu")){ 
      fertmax <- 1 } 
    
    ###################################################
    # Go through fert ids...
    for(ferti in 1:fertmax){
      ferti <<- ferti # make global
      print(paste(speciesName,"/",harvScen,"/",harvInten,"/ fert =",ferti))
      time0 <- Sys.time()
      outputAgei <-list()
      if(!(harvSceni%in%c("Powerline_under","Powerline_border"))) initAges <- c(NA,yearsToMem)
      initAgei <- init0 <- 1
      ##############################################################
      for(initAgei in init0:length(initAges)){ # 50#NA
        initAge <- initAges[initAgei]
        #if(initAgei==0){ 
        #  initAge <- NA
        #} else {initAge <- c(NA,yearsToMem)[initAgei]}
        toMem <- ls()
        dataS <- dataSorig
        if(CoeffSim){
          simInitData <- data.table(ba = 0.01, age = 0, dbh = 0.01, pine = 1, 
                                    spruce =  1, birch =  1, #decid =  species[4], 
                                    fert = ferti, h = 1.35, minpeat = minpeat0,
                                    landclass = landclass0, cons = 0)
          nnSim <- 1:nSitesRun0
          dataS$ba[nnSim] <- simInitData$ba
          dataS$age[nnSim] <- simInitData$age
          dataS$dbh[nnSim] <- simInitData$dbh
          if(speciesName%in%c("typical","kitu")){
            # Use different species compositions for all segments
            dataS$pine[nnSim] <- speciess[,"pine",speciesSeti,ferti]#simInitData$pine
            dataS$spruce[nnSim] <- speciess[,"spruce",speciesSeti,ferti]#simInitData$spruce
            dataS$birch[nnSim] <- 0*speciess[,"birch",speciesSeti,ferti]# simInitData$birch
            dataS$decid[nnSim] <- speciess[,"birch",speciesSeti,ferti]#simInitData$decid
          } else {
            # Use same average species compositions for all segments
            dataS$pine[nnSim] <- mean(speciess[,"pine",speciesSeti,ferti])#simInitData$pine
            dataS$spruce[nnSim] <- mean(speciess[,"spruce",speciesSeti,ferti])#simInitData$spruce
            dataS$birch[nnSim] <- 0# simInitData$birch
            dataS$decid[nnSim] <- mean(speciess[,"birch",speciesSeti,ferti])#simInitData$decid
          }
          
          if(speciesName=="typical"){
            dataS$fert[nnSim] <- ferttypical          
            #dataS$landclass[nnSim] <- simInitData$landclass
          } else if(speciesName=="kitu"){
            dataS$fert[nnSim] <- fertkitu          
            dataS$landclass[nnSim] <- 2
          } else {
            dataS$fert[nnSim] <- simInitData$fert
            dataS$landclass[nnSim] <- simInitData$landclass
          }
          
          dataS$h[nnSim] <- simInitData$h
          dataS$minpeat[nnSim] <- simInitData$minpeat
          dataS$cons[nnSim] <- simInitData$cons
          if(n0only) dataS <- dataS[nnSim,]
          #print(head(dataS))
          #ops <<- list(dataS)
        }
        
        #source_url("https://raw.githubusercontent.com/virpi-j/Hiilikartta/master/functions.R")
        #source("~/Hiilikartta/functions.R", local = T)
        if(manualRun){
          RCP=climScen; climdata=NULL; easyInit=FALSE; forceSaveInitSoil=F; cons10run = F; procDrPeat=F; outType = "hiiliKartta"; coeffPeat1=-240; coeffPeat2=70; coefCH4 = 0.34; coefN20_1 = 0.23; coefN20_2 = 0.077; landClassUnman=NULL; compHarvX = 0; funPreb = regionPrebas; initSoilCreStart=NULL; outModReStart=NULL; reStartYear=1; sampleX=dataS; P0currclim=NA; fT0=NA; sampleID <- 1
        }
        if(initAgei==1){
          out <- lapply(sampleIDs, function(jx) {
            runModel(jx,harvScen=harvScen, harvInten=harvInten, outType = "hiiliKartta", 
                     RCP = climScen, initAge = initAge, ingrowth = F, sampleX = dataS)})
          
          #if(harvSceni=="NoHarv"){ # for all harvest scenarios, the init state at age x comes from noharv-scenario
          if(is.na(initAge)){
            multiOut <- array(0,dim = c(dim(out[[1]]$restartMod$multiOut),length(sampleIDs)))
            GVOut <- array(0,dim = c(dim(out[[1]]$restartMod$GVout),length(sampleIDs)))
            reStartSoil <- array(0,dim = c(dim(out[[1]]$reStartSoil),length(sampleIDs)))
            for(ij in sampleIDs){ # for climate models
              if(length(dim(multiOut))==5){
                multiOut[,,,,ij] <- out[[ij]]$restartMod$multiOut
                GVOut[,,,ij] <- out[[ij]]$restartMod$GVout
                reStartSoil[,,,,ij] <- out[[ij]]$reStartSoil
              } else {
                multiOut[,,,,,ij] <- out[[ij]]$restartMod$multiOut
                GVOut[,,,ij] <- out[[ij]]$restartMod$GVout
                reStartSoil[,,,,,ij] <- out[[ij]]$reStartSoil
              }
              if(harvSceni=="NoHarv"){ # save clim data in the first run they are used
                climdata <- out[[ij]]$clim
                save(climdata, file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/HiiliKartta_climdata",r_no,"_clim",climModids[ij],".rdata"))
                print(paste("save climdata for",paste0("clim",climModids[ij])))
              }              
            }
            # average over climate models
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
            print(paste("Init stages saved for ages")); print(yearsToMem)
          }
          if(ingrowth){ # if ingrowth=T run again with ingrowth
            print("Run again with ingrowth!")
            out <- lapply(sampleIDs, function(jx) {
              print(paste("open data for clim",climModids[jx]))
              load(file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/HiiliKartta_climdata",r_no,"_clim",climModids[jx],".rdata"))
              runModel(jx,sampleX = dataS, harvScen=harvScen, harvInten=harvInten, outType = "hiiliKartta", 
                       ingrowth = ingrowth, climdata = climdata,
                       RCP = climScen, initAge = initAge)})
          }
          #}
        } else {
          out <- lapply(sampleIDs, function(jx) {
            print(paste("open data for clim",climModids[jx]))
            load(file=paste0("/scratch/project_2000994/PREBASruns/PREBAStesting/HiiliKartta_climdata",r_no,"_clim",climModids[jx],".rdata"))
            runModel(jx,sampleX = dataS, harvScen=harvScen, harvInten=harvInten, outType = "hiiliKartta", 
                     ingrowth = ingrowth,climdata = climdata, 
                     RCP = climScen, initAge = initAge)})
          #  tta <- aetgs
        }
        
        rm(list=setdiff(ls(),c(toMem,"out")))
        
        #if(initAgei>1){
        V <- Wround <- Wenergy <- litters <- ba <- H <- age <- nep <- wTot <- wGV <- soilC <- Vpine <- Vspruce <- Vbirch <- grossGrowth <-
          array(0,c(nSitesRun0,nYears,length(sampleIDs)),
                dimnames = list(paste0("site",1:nSitesRun0),
                                2014+1:nYears,
                                climMod[sampleIDs]))
        for(ij in 1:length(out)){
          V[,,ij] <- out[[ij]]$V
          Wround[,,ij] <- out[[ij]]$Wround
          Wenergy[,,ij] <- out[[ij]]$Wenergy
          ba[,,ij] <- out[[ij]]$ba
          H[,,ij] <- out[[ij]]$H
          age[,,ij] <- out[[ij]]$age
          nep[,,ij] <- out[[ij]]$nep
          wTot[,,ij] <- out[[ij]]$wTot
          wGV[,,ij] <- out[[ij]]$wGV
          soilC[,,ij] <- out[[ij]]$soilC
          litters[,,ij] <- out[[ij]]$litters
          Vpine[,,ij] <- out[[ij]]$Vpine
          Vspruce[,,ij] <- out[[ij]]$Vspruce
          Vbirch[,,ij] <- out[[ij]]$Vbirch
          grossGrowth[,,ij] <- out[[ij]]$grossGrowth
        }
        outputAgei[[initAgei]] <- list(V, Wround, Wenergy, ba, H, age, nep, wTot, wGV, soilC, litters, Vpine, Vspruce, Vbirch, grossGrowth)
        names(outputAgei[[initAgei]]) <- c("V","Wround", "Wenergy", "ba","H", "age", "nep", "wTot", "wGV", "soilC", "litters", "Vpine", "Vspruce", "Vbirch","grossGrowth")
        print(Sys.time()-time0)
        par(mfrow=c(3,1))
        plot(colMeans(V[,,1]),type="l", ylim = c(-0.5, max(V)),ylab="V", xlab="time", 
             main=paste(harvScen,"/",speciesName,"/ age",initAge,"/ V, fert",ferti))
        lines(colMeans(Vpine[,,1]),col="red")
        lines(colMeans(Vspruce[,,1]),col="blue")
        lines(colMeans(Vbirch[,,1]),col="green")
        plot(colMeans(H[,,1]),type="l", ylim = c(-0.5, max(H)),ylab="H", xlab="time",
             main=paste(harvScen,"/",speciesName,"/ age",initAge,"/ H, fert",ferti))
        plot(colMeans(grossGrowth[,,1]),type="l", 
             ylim = c(-0.5, max(grossGrowth)),ylab="grossgrowth", xlab="time",
             main=paste(harvScen,"/",speciesName,"/ age",initAge,"/ grossgrowth, fert",ferti))
        #}
      }
      output[[ferti]] <- outputAgei
      names(output[[ferti]]) <- paste0("initAge",inAs)
      print(paste0(harvScen," / fert",ferti," / age0 / V:"))
      print(output[[ferti]][[1]]$V[1,1:10,1])
      if(length(initAges)>1){
        print(paste0(harvScen," / fert",ferti," / age30 / V:"))
        print(output[[ferti]][[2]]$V[1,1:10,1])
      }
    }
    names(output) <- paste0("fert",1:fertmax)
    
    outFileePath <-"/scratch/project_2000994/PREBASruns/PREBAStesting/HiilikarttaResults/"
    #if(is.na(initAge)) initAge <- 0
    outFilee <- paste0(outFileePath,"HiiliKarttaTestPlots_rno",r_no,"_",harvScen,"_",harvInten,
                       "_species",speciesName,".pdf")
    if(toFile) pdf(outFilee)
    for(ij in 1:length(output[[1]][[1]])){
      par(mfrow=c(ceiling(sqrt(fertmax)),floor(sqrt(fertmax))))
      ymax <- max(output[[1]][[1]][[ij]])
      ymin <- min(output[[1]][[1]][[ij]])
      for(ferti in 1:fertmax){
        for(agei in 1:length(inAs)){
          ymax <- max(ymax,max(colMeans(output[[ferti]][[agei]][[ij]])))
          ymin <- min(ymin,min(colMeans(output[[ferti]][[agei]][[ij]])))
        }
      }
      ylims <- c(ymin,ymax)
      
      ageplot <- F
      if(ageplot){
        for(agei in 1:length(inAs)){
          par(mfrow=c(ceiling(sqrt(fertmax)),floor(sqrt(fertmax))))
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
                   main=paste("fert =",ferti,"init year",inAs[agei]))
              for(ik in sampleIDs){ # go through climate models
                #lines(as.numeric(names(Vmean)),Vi[,ij],col="blue")
                lines(ages[,ik],xi[,ik],col="blue")
              }
              if(ymin<0) lines(c(min(ages),max(ages)),c(0,0),col="black")
            }
          }
        }
      }
      varplot <- T
      if(toFile | varplot){
        for(agei in 1:length(inAs)){
          for(ferti in 1:fertmax){
            tmp <- output[[ferti]][[agei]][[ij]]
            xi <- apply(tmp,2:3,mean)
            xmean <- apply(tmp,2,mean)
            time <- 1:ncol(output[[ferti]][[agei]][[ij]])
            plot(time,xmean,xlab="time",ylab=names(output[[ferti]][[agei]])[ij],
                 type="l",ylim=ylims,lwd=2,
                 main=paste("fert =",ferti,"init year",inAs[agei]))
            for(ik in sampleIDs){ # go through climate models
              #lines(as.numeric(names(Vmean)),Vi[,ij],col="blue")
              lines(time,xi[,ik],col="blue")
            }
            if(ymin<0) lines(c(min(time),max(time)),c(0,0),col="black")
          }
        }
      }
    }
    if(toFile) dev.off()
    outFilee <- paste0(outFileePath,"HiiliKarttaTestPlots_rno",r_no,"_",harvScen,"_",harvInten,
                       "_species",speciesName,".rdata")
    save(output,file=outFilee)
    print(paste("Saved file",outFilee))
  }
}
###########
ij <- ij0 <- 1
for(ij in ij0:nsets){
  #species <<- speciess[ij,]
  speciesName <<- speciesNames[ij]
  runOut <- lapply(harvScens[1], function(jx) {
    runPerHarvScen(jx,speciesSeti = ij)})
  if(speciesName!="kitu"){ # for kitumaa, only noharv scen
    if(!parRuns){
      runOut <- lapply(harvScens[-1], function(jx) {
        runPerHarvScen(jx,speciesSeti = ij)})
    } else {
      runOut <- mclapply(harvScens[-1], function(jx) {
        runPerHarvScen(jx,speciesSeti = ij)
      }, mc.cores = 5, mc.silent=FALSE)      
    }      
  }
}


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

Sys.chmod(list.dirs("initSoilC"), "0777",use_umask=FALSE)
f <- list.files("initSoilCunc", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)

Sys.chmod(list.dirs("plots"), "0777",use_umask=FALSE)
f <- list.files("plots", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)
