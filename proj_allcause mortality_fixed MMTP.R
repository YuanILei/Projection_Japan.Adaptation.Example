################################################################################
# Aim: to summary & plot country- and prefecture-specific projected 
# temperature data obtained and processed by Paul
# 3rd_stage_proj, all-cause mortality
# ADAPTATION SENARIO: CONSTANT MMTP
# Created and Updated by Lei Yuan on April 20th, 2024
################################################################################
rm(list=ls())
setwd("~/Desktop/data/Japan/ISIMIP3b/")
library(splines)
library(dlnm)
library(mixmeta)
library(MASS)
library(lubridate)
################################################################################
# IMPORT DATA
# BLUP values
load("blup_allmorta_1519.rda") # baseline ERF coef and vcov

# calibrated modeled temperatures 
load("isimip3b_tasAdjust_4rcp_list.rda")# nlist--calibrated modeled temperatures
load("jpmortality_7219_list.rda") # dlist--historical temperature

# csv file--prefecture-specific mmt,mmtp with blup (variable: "prefecture","mmtp","mmt")
pref <- read.csv("mintemppref_allmorta.csv",stringsAsFactors = FALSE)
################################################################################
# create objects
tseq1 <- seq(as.Date("2015-01-01"),as.Date("2019-12-31"),"day") #for historical data
tseq2 <- seq(as.Date("2010-01-01"),as.Date("2099-12-31"),"day") #for projected temperatures
# remove Feb 29th
tseq1 <- tseq1[-grep("02-29",tseq1)]
tseq2 <- tseq2[-grep("02-29",tseq2)]

gcm <- c("gfdl.esm4","ipsl.cm6a.lr","mpi.esm1.2.hr","mri.esm2.0","ukesm1.0.ll") # 5 GCMs
ssp <- c("ssp126","ssp245","ssp370","ssp585") # 4 SSPs
period <- paste(201:209,0,"-",substr(201:209,3,3),9,sep="") # by decades
nsim <- 1000

# model specs
varfun <- "ns"
varper <- c(10,75,90)
# vardf <- 4 # ns+3 knots

# create arrays for data storage
# attributable numbers and fractions by prefecture
anpref <- afpref <- array(0,dim=c(nrow(pref),3,2,3,length(period),length(gcm)+1,length(ssp)),
                          dimnames=list(pref$prefecture,c("est","ci.l","ci.u"),
                                        c("abs","rel"),c("tot","cold","heat"),
                                        period,c("ensemble",gcm),ssp)) 
dimnames(anpref)
# attributable numbers and fractions for entire Japan
anctry <- afctry <- array(0,dim=c(3,2,3,length(period),length(gcm)+1,length(ssp)),
                          dimnames=list(c("est","ci.l","ci.u"),
                                        c("abs","rel"),c("tot","cold","heat"),
                                        period,c("ensemble",gcm),ssp))
dimnames(anctry)
# temperature
tyr <- array(NA,c(nrow(pref),length(2010:2099),length(gcm),length(ssp)),
             dimnames=list(pref$prefecture,2010:2099,gcm,ssp))
dimnames(tyr)
# save centering of the temperatures by decade
cen.df <- data.frame(matrix(NA,nrow=nrow(pref),ncol=length(period)+2)) 
colnames(cen.df) <- c("pref.name","pref.code",period)
cen.df[,c("pref.name","pref.code")] <- pref[,c("prefname","prefecture")]
# vector for projected outcome
dperiod <- rep(0,nrow(pref))
names(dperiod) <- pref$prefname
################################################################################
# loop for mortality
# i=j=k=p=1L #for testing
for (i in seq(ssp)) {
  # print
  cat("\n\n",ssp[i],"\n")
  # STORE UNCERTAINTY
  # DIMENSION - ABSOLUTE AN/DIFFERENCE IN AN
  anprefsim <- array(0,dim=c(nrow(pref),2,3,length(period),length(gcm),nsim+1),
                     dimnames=list(pref$prefecture,c("abs","rel"),c("tot","cold","heat"),
                                   period,gcm,c("est",paste0("sim",seq(nsim)))))
  
  # loop by prefecture
  for (j in pref$prefecture){
    # print
    cat(j," ")
    # tmean projection
    tproj1 <- nlist[[j]] # temp projection by prefecture
    tproj2 <- tproj1[tproj1$date%in%tseq2,c(1,grep(ssp[i],colnames(tproj1)))] #grep ssp
    tyr[j,,,i] <- apply(tproj2[,-1],2,tapply,substr(tproj2$date,1,4),mean) #take mean by year
    # outcome projection
    sdat <- dlist[[j]]
    sdat <- sdat[sdat$date[-grep("-02-29",sdat$date)],]
    sdat <- sdat[sdat$date%in%tseq1,]  #Feb 29 is removed
    # -----------------------------------------------
    # 1---baseline sub-period outcome counts
    # all pop
    ddoy <- tapply(sdat$all,as.numeric(format(sdat$date,"%j")),mean)[seq(365)] #average daily counts over baseline subperiods
    # 
    dproj <- rep(ddoy,length=length(tseq2))
    speriod <- factor(rep(period,each=365*10))
    dperiod[j] <- sum(ddoy)*10
    # argvar and centre
    argvar <- list(fun=varfun,knots=quantile(sdat$tmean,varper/100,na.rm=T),
                   Bound=range(sdat$tmean,na.rm=T))
    # parameters
    coef <- blup[[j]]$blup
    vcov <- blup[[j]]$vcov
    # multivariate normal distribution
    set.seed(1515) #change seed
    coefsim <- mvrnorm(nsim,coef,vcov)
    
    # loop by GCM
    for (k in seq(gcm)) {
      # scenario name
      sce <- paste0(gcm[k],"_",ssp[i])
      
      # loop by period / decade
      for (p in seq(period)) {
        yrs <- as.integer(substr(period[p],1,4)):as.integer(paste0("20",substr(period[p],6,8))) # years in decade
        # get centered basis # 2--full adaptation: fixed MMTP
        bvar <- do.call(onebasis,c(list(x=tproj2[year(tproj2$date)%in%yrs,sce]),argvar))
        cen <- quantile(tproj2[year(tproj2$date)%in%yrs,sce],pref$mmtp[j]/100,na.rm=T)
        cen.df[j,period[p]] <- cen # store temperature
        cenvec <- do.call(onebasis,c(list(x=cen),argvar))
        bvarcen <- scale(bvar,center=cenvec,scale=F)
        # indicator for heat
        indheat <- tproj2[year(tproj2$date)%in%yrs,sce]>cen
        # daily attributable outcomes
        an <- (1-exp(-bvarcen%*%coef))*dproj[which(year(tproj2$date)%in%yrs)]
        # sum by range and period, getting mean
        anprefsim[j,1,1,p,k,1] <- sum(an) #both cold and heat
        anprefsim[j,1,2,p,k,1] <- sum(an[!indheat]) #cold
        anprefsim[j,1,3,p,k,1] <- sum(an[indheat]) #heat
        
        # loop by simulation
        for (s in seq(nsim)){
          # daily attributable outcomes
          an <- (1-exp(-bvarcen%*%coefsim[s,]))*dproj[which(year(tproj2$date)%in%yrs)]
          # COMPUTE THE RESULTS FOR EACH RANGE AND PERIOD AND SUM
          # NB: ACCOUNT FOR NO TEMPERATURE BELOW/ABOVE CEN FOR GIVEN PERIODS
          anprefsim[j,1,1,p,k,s+1] <- sum(an)
          anprefsim[j,1,2,p,k,s+1] <- sum(an[!indheat])
          anprefsim[j,1,3,p,k,s+1] <- sum(an[indheat])
        }
      }
    }
  }
  
  
  ################################################################################
  # RELATIVE TO 2010-19
  anprefsim[,2,,,,] <- anprefsim[,1,,,,] - anprefsim[,1,,rep(1,9),,]
  
  ################################################################################
  # AGGREGATE CITY-SPECIFIC DATA AT HIGHER LEVEL
  # aggregate to country-level
  anctrysim <- apply(anprefsim[,,,,,],2:6,sum)
  
  ################################################################################
  # ATTRIBUTABLE NUMBERS, BY GCM
  # prefecture-level, WITH CI
  anpref[,1,,,,-1,i] <- apply(anprefsim[,,,,,1],1:5,mean)
  anpref[,2,,,,-1,i] <- apply(anprefsim[,,,,,-1],1:5,quantile,0.025)
  anpref[,3,,,,-1,i] <- apply(anprefsim[,,,,,-1],1:5,quantile,0.975)
  
  # country level, WITH CI
  anctry[1,,,,-1,i] <- apply(anctrysim[,,,,1],1:4,mean)
  anctry[2,,,,-1,i] <- apply(anctrysim[,,,,-1],1:4,quantile,0.025)
  anctry[3,,,,-1,i] <- apply(anctrysim[,,,,-1],1:4,quantile,0.975)
  
  ################################################################################
  # ATTRIBUTABLE NUMBERS, ENSEMBLE
  # prefecture-level, WITH CI
  anpref[,1,,,,1,i] <- apply(anprefsim[,,,,,1],1:4,mean)
  anpref[,2,,,,1,i] <- apply(anprefsim[,,,,,-1],1:4,quantile,0.025)
  anpref[,3,,,,1,i] <- apply(anprefsim[,,,,,-1],1:4,quantile,0.975)
  
  # country level, WITH CI
  anctry[1,,,,1,i] <- apply(anctrysim[,,,,1],1:3,mean)
  anctry[2,,,,1,i] <- apply(anctrysim[,,,,-1],1:3,quantile,0.025)
  anctry[3,,,,1,i] <- apply(anctrysim[,,,,-1],1:3,quantile,0.975) 
  
  ################################################################################
  # ATTRIBUTABLE FRACTIONS (ALSO COMPUTING THE NET EFFECT)
  afpref[,,,,,,i] <- anpref[,,,,,,i]/dperiod*100
  afctry[,,,,,i] <- anctry[,,,,,i]/sum(dperiod)*100
}

################################################################################
# CHANGE IN AF (%) BETWEEN 2010-19 & 2090-99 (RCP8.5)

# HEAT
afctry[,,"heat","2090-99","ensemble","ssp585"]

# COLD
afctry[,,"cold","2090-99","ensemble","ssp585"]

anpref
afpref
afctry
anctry
tyr
################################################################################

# SAVE 
# all pop
save(list=c("cen.df","anpref","afpref","anctry","afctry","tyr"),file="temp_attr_allmort_nopop_fixed mmtp_all.rda")



