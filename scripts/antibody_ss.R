# Sample size calculation for Vacc-iNTS
## Antibody as endpoint

library(ggplot2)
library(RColorBrewer)
#library(lme4)
#library(lmerTest)
library(gee) # using GEE instead of LMM; fixed intercept in LMM captured some of the effect due to exposure...
library(doParallel)
library(pwr)
library(here)
library(rio)
library(Hmisc)

args<-commandArgs(T)

if(!is.na(args[1])){
  cat(paste(sep="","Using ",args[1]," cores!\n\n"))
  registerDoParallel(cores=args[1])
}else{
  cat(paste(sep="","Using ",2," cores!\n\n"))
  registerDoParallel(cores=2)
}

set.seed(20200128)

# power parameters
alpha<-0.05
power<-0.8
B<-1e3 # number of simulations

## Malaria infection

if(is.na(args[2]) | args[2]=="malaria"){
  # simulate data
  simulateDataMalAb<-function(n,prevBin=0.16,abSD=30,abSlopeMonth=2.5,iccSD=0.03,effectSize=-12.5,effectSizeSD=30,pDropOut=0.2,ageMin=8,ageMax=36){
    # n<-100 # this would yield 2n measurements, each pair of measurements for the n individuals 3 months apart
    # prevBin<-0.2 # malaria prevalence in children 6-60 months, Southern Malawi, MIS 2017 p.43 Fig 4.4, https://dhsprogram.com/pubs/pdf/MIS28/MIS28.pdf
    # abSD<- 39.8 # Nyirenda et al (2014)
    # abSlopeMonth<-(3.2) # Nyirenda et al (2014)
    # ageMin<-8; ageMax<-28 # Nyirenda et al (2014)
    
    pid<-paste(sep="","p",1:n)
    ages<-runif(n=n,min=ageMin,max=ageMax) 
    binCov<-rbinom(n=2*n,size=1,prob=prevBin) # allows participants to change state for each visit; assumed independent risk at each visit and risk the same at all ages -- both of which may not be entirely realistic...
    
    u<-rnorm(n=n,mean=0,sd=iccSD)
    v<-rnorm(n=n,mean=0,sd=iccSD)
    ab8m<-rnorm(n=n,mean=0,sd=abSD)
    
    datSim<-data.frame(pid=rep(pid,2),
                       visit=c(rep(1,n),rep(2,n)),
                       age=c(ages,ages+rnorm(n,mean=3,sd=0.05)), # approximately 3 months between initial recruitment and follow-up but allowing minor deviations from planned visit dates
                       ageSlopeEffectLimit=NA,
                       exposure=ifelse(binCov==1,"exposed","unexposed"),
                       u=rep(u,2),
                       v=rep(v,2),
                       ab8m=rep(ab8m,2),
                       ab=NA)
    
    datSim$exposure<-factor(datSim$exposure,levels=c("unexposed","exposed"))
    datSim$agem8<-datSim$age-8
    
    for(i in 1:nrow(datSim)){
      pid<-datSim$pid[i]
      idx1<-which(datSim$pid==pid & datSim$visit==1)
      idx2<-which(datSim$pid==pid & datSim$visit==2)
      datSim$ageSlopeEffectLimit[idx2]<-ifelse(datSim$exposure[idx1]=="exposed",datSim$age[idx2]-0.5,datSim$age[idx2])
      #datSim$ageSlopeEffectLimit[idx1]<-datSim$age[idx1]
      
    }
    
    datSim$ab<-datSim$ab8m+datSim$u+(abSlopeMonth+datSim$v)*(datSim$ageSlopeEffectLimit-8)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
    datSim$ab[datSim$visit==2 & rbinom(n=n,size=1,prob=pDropOut)==1]<-NA
    datSim$ab[datSim$ab<=0]<-0.0001
    
    return(datSim)
  }
  
  # detect a change in slope of ab with age for i) malaria cases, ii) malaria vaccinated, iii) other binary covariates (linear mixed model approach)
  fitGEETestCoeff<-function(dat,alpha=0.05,coefName1="exposureexposed",coefName2=NULL){
    library(gee)
    mod<-gee(ab~agem8*exposure,id=as.factor(pid),data=dat[order(dat$pid,dat$visit),],corstr="exchangeable",na.action=na.omit)
    p1<-2*(1-pnorm(abs(summary(mod)$coefficients[coefName1,"Robust z"])))
    p2<-2*(1-pnorm(abs(summary(mod)$coefficients[coefName2,"Robust z"])))
    
    if(p1<alpha){
      res<-1
    }else{
      res<-0
    }
    
    if(!is.null(coefName2)){
      if(p2<alpha){
        res<-c(res,1)
      }else{
        res<-c(res,0)
      }
    }
    
    return(res)
  }
}
  
  
  # malaria infection; prevalence ~16%; try sample sizes 2000 and effect sizes 0.1, 0.2, 0.3
  gr<-expand.grid(c(2000, 2500, 3000),c(6,9,12,15,18),c(0.2))
  #gr<-expand.grid(c(2500),c(12),c(0.2))
  datPow<-data.frame(
    alpha=0.05,
    malariaPrev=0.16,
    effectSize=gr[,2],
    n=gr[,1],
    pDropOut=gr[,3],
    powerShift=NA,
    powerSlope=NA
  )
  
  for(j in 1:nrow(datPow)){
    cat(paste(sep="","## effectSize = ",datPow$effectSize[j],", n = ",datPow$n[j],", pDropOut = ",datPow$pDropOut[j]," ##\n"))
    
  
    res<-foreach(i = 1:B) %dopar% {
     
      datSim<-simulateDataMalAb(n=datPow$n[j],prevBin=datPow$malariaPrev[j],abSD=30,abSlopeMonth=(2.5),iccSD=0.03,effectSize=datPow$effectSize[j],effectSizeSD=30,pDropOut=datPow$pDropOut[j],ageMin=8,ageMax=36)
      tmp<-fitGEETestCoeff(dat=datSim,coefName1="exposureexposed",coefName2="agem8:exposureexposed",alpha=datPow$alpha[j])
     
      tmp
    }
    
    tmp<-matrix(unlist(res),byrow=T,ncol=2)
    resInt<-tmp[,1]
    resSlp<-tmp[,2]
    
    datPow$powerShift[j]<-sum(resInt)/B
    datPow$powerSlope[j]<-sum(resSlp)/B
    cat(paste(sep="","power to detect a shift in ab = ",datPow$powerShift[j],", power to detect a change in ab-age slope = ",datPow$powerSlope[j]," [",datPow$effectSize[j],", ",datPow$n[j],", ",datPow$pDropOut[j],"]\n\n"))
  }
  
datPowMalInfAb2 <- datPow
export(datPow, here("data", "DawPowMalInf_Ab2.RData"))
  
  # produce power curves
pdf(file="~/work/PhD/VacciNTS/Sample size/Vacc-iNTS_SampleSize_MalariaInfection_powerCurves_pDO0.2.pdf",width=16,height=9)
  ggplot(data=datPowMalInfAb2[datPowMalInfAb2$pDropOut==0.2,],mapping=aes(x=effectSize,y=powerShift,group=as.factor(n),color=as.factor(n))) +
    geom_point(size=3) +
    geom_line(lwd=1.5) +
    geom_hline(yintercept=0.8,lwd=1,lty=2,col="darkgrey") + 
    scale_color_manual(values=colorRampPalette(brewer.pal(8,"Spectral"))(length(levels(factor(datPow$n)))),name="Number of paired samples") +
    scale_x_continuous(breaks = seq(0,20,1)) +
    ylab("Power") +
    xlab("Effect size (shift in average ab)") +
    ggtitle("Power curves: Impact of malaria on NTS antibody responses, drop-out (20%)") +
    theme(text=element_text(size=16)) +
    theme_bw()
  dev.off()
  
  # produce a few example figures of data generated
  c1<-col2rgb("steelblue")
  c1<-rgb(c1[1],c1[2],c1[3],alpha=75,maxColorValue=255)
  c2<-col2rgb("orange")
  c2<-rgb(c2[1],c2[2],c2[3],alpha=75,maxColorValue=255)
  
  pdf(file="Vacc-iNTS_SampleSize_MalariaInfection_n1000_es0.1_pDO0.2.pdf",width=16,height=9)
  datSim<-simulateDataMalAb(n=1000,prevBin=0.16,abSD=30,abSlopeMonth=2.5,iccSD=0.03,effectSize=-12.5,effectSizeSD=30,pDropOut=0.2,ageMin=8,ageMax=36)
  mod<-gee(ab~agem8*exposure,id=as.factor(pid),data=datSim[order(datSim$pid,datSim$visit),],corstr="exchangeable",na.action=na.omit)
  datSim$pid_n <- as.numeric(stringr::str_extract(datSim$pid, "\\d+"))
  mod<-gee(ab~agem8*exposure,id=pid_n,data=datSim,corstr="exchangeable")
  
  g<-ggplot(data=datSim,mapping=aes(x=agem8,y=ab,col=exposure, group = pid_n)) +
    geom_point() + 
    geom_line() +
    scale_color_manual(values=c(c1,c2)) +
    geom_abline(lwd=1,intercept=coef(mod)["(Intercept)"]+coef(mod)["exposureexposed"],slope=coef(mod)["agem8"]+coef(mod)["agem8:exposureexposed"],col="orange") + 
    geom_abline(lwd=1,intercept=coef(mod)["(Intercept)"],slope=coef(mod)["agem8"],col="steelblue") +
    theme_bw() +
    xlab("Age in months") +
    ylab("Antibody titre") +
    ggtitle("Example results for one simulation: antibody responses with age by malaria exposure")
  print(g)
  dev.off()

  
 ############################
  ## Anaemia
  
  
  if(is.na(args[2]) | args[2]=="anaemia"){
    # simulate data
    simulateDataAnAb<-function(n,prevBin=0.2,abSD=30,abSlopeMonth=2.5,iccSD=0.03,effectSize=-12.5,effectSizeSD=30,pDropOut=0.2,ageMin=8,ageMax=36){
      # n<-100 # this would yield 2n measurements, each pair of measurements for the n individuals 3 months apart
      # prevBin<-0.2 # malaria prevalence in children 6-60 months, Southern Malawi, MIS 2017 p.43 Fig 4.4, https://dhsprogram.com/pubs/pdf/MIS28/MIS28.pdf
      # abSD<- 39.8 # Nyirenda et al (2014)
      # abSlopeMonth<-(3.2) # Nyirenda et al (2014)
      # ageMin<-8; ageMax<-28 # Nyirenda et al (2014)
      
      pid<-paste(sep="","p",1:n)
      ages<-runif(n=n,min=ageMin,max=ageMax) 
      binCov<-rbinom(n=2*n,size=1,prob=prevBin) # allows participants to change state for each visit; assumed independent risk at each visit and risk the same at all ages -- both of which may not be entirely realistic...
      
      u<-rnorm(n=n,mean=0,sd=iccSD)
      v<-rnorm(n=n,mean=0,sd=iccSD)
      ab8m<-rnorm(n=n,mean=0,sd=abSD)
      
      datSim<-data.frame(pid=rep(pid,2),
                         visit=c(rep(1,n),rep(2,n)),
                         age=c(ages,ages+rnorm(n,mean=3,sd=0.05)), # approximately 3 months between initial recruitment and follow-up but allowing minor deviations from planned visit dates
                         ageSlopeEffectLimit=NA,
                         exposure=ifelse(binCov==1,"exposed","unexposed"),
                         u=rep(u,2),
                         v=rep(v,2),
                         ab8m=rep(ab8m,2),
                         ab=NA)
      
      datSim$exposure<-factor(datSim$exposure,levels=c("unexposed","exposed"))
      datSim$agem8<-datSim$age-8
      
      for(i in 1:nrow(datSim)){
        pid<-datSim$pid[i]
        idx1<-which(datSim$pid==pid & datSim$visit==1)
        idx2<-which(datSim$pid==pid & datSim$visit==2)
        #datSim$ageSlopeEffectLimit[i]<-ifelse(datSim$exposure[idx1]=="exposed",datSim$age[idx1]-0.5,datSim$age[i])
        datSim$ageSlopeEffectLimit[idx2]<-ifelse(datSim$exposure[idx1]=="exposed",datSim$age[idx2]-0.5,datSim$age[idx2])
      }
      
      #datSim$ab<-datSim$ab8m+datSim$u+(abSlopeMonth+datSim$v)*(datSim$age-8-ifelse(datSim$visit==1 & datSim$exposure=="exposed",1,ifelse(datSim$visit==2 & datSim$exposure=="exposed",0.5,0)))+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
      datSim$ab<-datSim$ab8m+datSim$u+(abSlopeMonth+datSim$v)*(datSim$ageSlopeEffectLimit-8)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
      #datSim$ab<-datSim$ab8m+datSim$u+(abSlopeMonth+datSim$v)*(datSim$age-8)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
      datSim$ab[datSim$visit==2 & rbinom(n=n,size=1,prob=pDropOut)==1]<-NA
      datSim$ab[datSim$ab<=0]<-0.0001
      
      return(datSim)
    }
    
    # detect a change in slope of ab with age for i) malaria cases, ii) malaria vaccinated, iii) other binary covariates (linear mixed model approach)
    fitGEETestCoeff<-function(dat,alpha=0.05,coefName1="exposureexposed",coefName2=NULL){
      library(gee)
      #mod<-lmer(ab~agem8*exposure+(1|pid),REML=F,data=dat) # REML=F as we test fixed parameters; cannot have a second random effect even though this is how data was generated; would need more data points per participant...
      mod<-gee(ab~agem8*exposure,id=as.factor(pid),data=dat[order(dat$pid,dat$visit),],corstr="exchangeable",na.action=na.omit)
      p1<-2*(1-pnorm(abs(summary(mod)$coefficients[coefName1,"Robust z"])))
      p2<-2*(1-pnorm(abs(summary(mod)$coefficients[coefName2,"Robust z"])))
      
      #if(summary(mod)$coefficients[coefName1,"Pr(>|t|)"]<alpha){
      if(p1<alpha){
        res<-1
      }else{
        res<-0
      }
      
      if(!is.null(coefName2)){
        #if(summary(mod)$coefficients[coefName2,"Pr(>|t|)"]<alpha){
        if(p2<alpha){
          res<-c(res,1)
        }else{
          res<-c(res,0)
        }
      }
      
      return(res)
    }
  }
  
  
  # Anaemia; prevalence ~20%; try sample sizes 2000 and effect sizes 0.1, 0.2, 0.3
  # gr<-expand.grid(c(2500),c(12),c(0.2))
  gr<-expand.grid(c(2000, 2500, 3000),c(6,9,12,15,18),c(0.2))
  datPow<-data.frame(
    alpha=0.05,
    AnaemPrev=0.2,
    effectSize=gr[,2],
    n=gr[,1],
    pDropOut=gr[,3],
    powerShift=NA,
    powerSlope=NA
  )
  
  for(j in 1:nrow(datPow)){
    cat(paste(sep="","## effectSize = ",datPow$effectSize[j],", n = ",datPow$n[j],", pDropOut = ",datPow$pDropOut[j]," ##\n"))
    
    #resInt<-rep(NA,B)
    #resSlp<-rep(NA,B)
    
    res<-foreach(i = 1:B) %dopar% {
      #if(i %% 100 == 0){
      #	cat(paste(sep="","..",i,"\n"))
      #}
      
      datSim<-simulateDataAnAb(n=datPow$n[j],prevBin=datPow$AnaemPrev[j],abSD=30,abSlopeMonth=(2.5),iccSD=0.03,effectSize=datPow$effectSize[j],effectSizeSD=30,pDropOut=datPow$pDropOut[j],ageMin=8,ageMax=36)
      tmp<-fitGEETestCoeff(dat=datSim,coefName1="exposureexposed",coefName2="agem8:exposureexposed",alpha=datPow$alpha[j])
      #resInt[i]<-tmp[1]
      #resSlp[i]<-tmp[2]
      tmp
    }
    
    tmp<-matrix(unlist(res),byrow=T,ncol=2)
    resInt<-tmp[,1]
    resSlp<-tmp[,2]
    
    datPow$powerShift[j]<-sum(resInt)/B
    datPow$powerSlope[j]<-sum(resSlp)/B
    cat(paste(sep="","power to detect a shift in ab = ",datPow$powerShift[j],", power to detect a change in ab-age slope = ",datPow$powerSlope[j]," [",datPow$effectSize[j],", ",datPow$n[j],", ",datPow$pDropOut[j],"]\n\n"))
  }
  
  # save(list=ls(),file="Vacc-iNTS_SampleSize_MalariaInfection.RData")
  datPowAnaemAb <- datPow
  export(datPow, here("data", "datPowAnaemAb.RData"))
  
  # produce power curves
  pdf(file="~/work/PhD/VacciNTS/Sample size/Vacc-iNTS_SampleSize_MalariaInfection_powerCurves_pDO0.2.pdf",width=16,height=9)
  ggplot(data=datPowAnaemAb[datPowAnaemAb$pDropOut==0.2,],mapping=aes(x=effectSize,y=powerShift,group=as.factor(n),color=as.factor(n))) +
    geom_point(size=3) +
    geom_line(lwd=1.5) +
    geom_hline(yintercept=0.8,lwd=1,lty=2,col="darkgrey") + 
    scale_color_manual(values=colorRampPalette(brewer.pal(8,"Spectral"))(length(levels(factor(datPow$n)))),name="Number of paired samples") +
    scale_x_continuous(breaks = seq(0,20,1)) +
    ylab("Power") +
    xlab("Effect size (shift in average ab)") +
    ggtitle("Power curves: impact of anaemia on NTS antibody responses, assuming 20% drop-out.") +
    theme(text=element_text(size=16)) +
    theme_bw()
  dev.off()
  
  # produce a few example figures of data generated
  c1<-col2rgb("steelblue")
  c1<-rgb(c1[1],c1[2],c1[3],alpha=75,maxColorValue=255)
  c2<-col2rgb("orange")
  c2<-rgb(c2[1],c2[2],c2[3],alpha=75,maxColorValue=255)
  
  pdf(file="Vacc-iNTS_SampleSize_MalariaInfection_n1000_es0.1_pDO0.2.pdf",width=16,height=9)
  datSim<-simulateDataAnAb(n=1000,prevBin=0.2,abSD=30,abSlopeMonth=2.5,iccSD=0.03,effectSize=-12.5,effectSizeSD=30,pDropOut=0.2,ageMin=8,ageMax=36)
  datSim$pid_n <- as.numeric(stringr::str_extract(datSim$pid, "\\d+"))
  mod<-gee(ab~agem8*exposure,id=pid_n,data=datSim,corstr="exchangeable")
  
  g<-ggplot(data=datSim,mapping=aes(x=agem8,y=ab,col=exposure, group = pid_n)) +
    geom_point() + 
    geom_line() +
    scale_color_manual(values=c(c1,c2)) +
    geom_abline(lwd=1,intercept=coef(mod)["(Intercept)"]+coef(mod)["exposureexposed"],slope=coef(mod)["agem8"]+coef(mod)["agem8:exposureexposed"],col="orange") + 
    geom_abline(lwd=1,intercept=coef(mod)["(Intercept)"],slope=coef(mod)["agem8"],col="steelblue") +
    theme_bw() +
    xlab("Age in months") +
    ylab("Antibody titre") +
    ggtitle("Example results for one simulation: antibody responses with age by malaria exposure")
  print(g)
  dev.off()

###############################################
## Malnutrition
  
  if(is.na(args[2]) | args[2]=="malnutrition"){
    # simulate data
    simulateDataMalnAb<-function(n,prevBin=0.02,abSD=30,abSlopeMonth=2.5,iccSD=0.03,effectSize=-12.5,effectSizeSD=30,pDropOut=0.2,ageMin=8,ageMax=36){
      # n<-100 # this would yield 2n measurements, each pair of measurements for the n individuals 3 months apart
      # prevBin<-0.2 # malaria prevalence in children 6-60 months, Southern Malawi, MIS 2017 p.43 Fig 4.4, https://dhsprogram.com/pubs/pdf/MIS28/MIS28.pdf
      # abSD<- 39.8 # Nyirenda et al (2014)
      # abSlopeMonth<-(3.2) # Nyirenda et al (2014)
      # ageMin<-8; ageMax<-28 # Nyirenda et al (2014)
      
      pid<-paste(sep="","p",1:n)
      ages<-runif(n=n,min=ageMin,max=ageMax) 
      binCov<-rbinom(n=2*n,size=1,prob=prevBin) # allows participants to change state for each visit; assumed independent risk at each visit and risk the same at all ages -- both of which may not be entirely realistic...
      
      u<-rnorm(n=n,mean=0,sd=iccSD)
      v<-rnorm(n=n,mean=0,sd=iccSD)
      ab8m<-rnorm(n=n,mean=0,sd=abSD)
      
      datSim<-data.frame(pid=rep(pid,2),
                         visit=c(rep(1,n),rep(2,n)),
                         age=c(ages,ages+rnorm(n,mean=3,sd=0.05)), # approximately 3 months between initial recruitment and follow-up but allowing minor deviations from planned visit dates
                         ageSlopeEffectLimit=NA,
                         exposure=ifelse(binCov==1,"exposed","unexposed"),
                         u=rep(u,2),
                         v=rep(v,2),
                         ab8m=rep(ab8m,2),
                         ab=NA)
      
      datSim$exposure<-factor(datSim$exposure,levels=c("unexposed","exposed"))
      datSim$agem8<-datSim$age-8
      
      for(i in 1:nrow(datSim)){
        pid<-datSim$pid[i]
        idx1<-which(datSim$pid==pid & datSim$visit==1)
        idx2<-which(datSim$pid==pid & datSim$visit==2)
        #datSim$ageSlopeEffectLimit[i]<-ifelse(datSim$exposure[idx1]=="exposed",datSim$age[idx1]-0.5,datSim$age[i])
        datSim$ageSlopeEffectLimit[idx2]<-ifelse(datSim$exposure[idx1]=="exposed",datSim$age[idx2]-0.5,datSim$age[idx2])
      }
      
      #datSim$ab<-datSim$ab8m+datSim$u+(abSlopeMonth+datSim$v)*(datSim$age-8-ifelse(datSim$visit==1 & datSim$exposure=="exposed",1,ifelse(datSim$visit==2 & datSim$exposure=="exposed",0.5,0)))+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
      datSim$ab<-datSim$ab8m+datSim$u+(abSlopeMonth+datSim$v)*(datSim$ageSlopeEffectLimit-8)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
      #datSim$ab<-datSim$ab8m+datSim$u+(abSlopeMonth+datSim$v)*(datSim$age-8)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
      datSim$ab[datSim$visit==2 & rbinom(n=n,size=1,prob=pDropOut)==1]<-NA
      datSim$ab[datSim$ab<=0]<-0.0001
      
      return(datSim)
    }
    
    # detect a change in slope of ab with age for i) malaria cases, ii) malaria vaccinated, iii) other binary covariates (linear mixed model approach)
    fitGEETestCoeff<-function(dat,alpha=0.05,coefName1="exposureexposed",coefName2=NULL){
      library(gee)
      #mod<-lmer(ab~agem8*exposure+(1|pid),REML=F,data=dat) # REML=F as we test fixed parameters; cannot have a second random effect even though this is how data was generated; would need more data points per participant...
      mod<-gee(ab~agem8*exposure,id=as.factor(pid),data=dat[order(dat$pid,dat$visit),],corstr="exchangeable",na.action=na.omit)
      p1<-2*(1-pnorm(abs(summary(mod)$coefficients[coefName1,"Robust z"])))
      p2<-2*(1-pnorm(abs(summary(mod)$coefficients[coefName2,"Robust z"])))
      
      #if(summary(mod)$coefficients[coefName1,"Pr(>|t|)"]<alpha){
      if(p1<alpha){
        res<-1
      }else{
        res<-0
      }
      
      if(!is.null(coefName2)){
        #if(summary(mod)$coefficients[coefName2,"Pr(>|t|)"]<alpha){
        if(p2<alpha){
          res<-c(res,1)
        }else{
          res<-c(res,0)
        }
      }
      
      return(res)
    }
  }
  
  
  # malnutritionprevalence ~2%; try sample sizes 2000 and effect sizes 0.1, 0.2, 0.3
  # gr<-expand.grid(c(2500),c(12),c(0.2))
  gr<-expand.grid(c(2500),c(40,50),c(0.2))
  #gr<-expand.grid(c(2000, 2500, 3000),c(30,40,50),c(0.2))
  datPow<-data.frame(
    alpha=0.05,
    malnPrev=0.02,
    effectSize=gr[,2],
    n=gr[,1],
    pDropOut=gr[,3],
    powerShift=NA,
    powerSlope=NA
  )
  
  for(j in 1:nrow(datPow)){
    cat(paste(sep="","## effectSize = ",datPow$effectSize[j],", n = ",datPow$n[j],", pDropOut = ",datPow$pDropOut[j]," ##\n"))
    
    #resInt<-rep(NA,B)
    #resSlp<-rep(NA,B)
    
    res<-foreach(i = 1:B) %dopar% {
      #if(i %% 100 == 0){
      #	cat(paste(sep="","..",i,"\n"))
      #}
      
      datSim<-simulateDataMalnAb(n=datPow$n[j],prevBin=datPow$malnPrev[j],abSD=30,abSlopeMonth=(2.5),iccSD=0.03,effectSize=datPow$effectSize[j],effectSizeSD=30,pDropOut=datPow$pDropOut[j],ageMin=8,ageMax=36)
      tmp<-fitGEETestCoeff(dat=datSim,coefName1="exposureexposed",coefName2="agem8:exposureexposed",alpha=datPow$alpha[j])
      #resInt[i]<-tmp[1]
      #resSlp[i]<-tmp[2]
      tmp
    }
    
    tmp<-matrix(unlist(res),byrow=T,ncol=2)
    resInt<-tmp[,1]
    resSlp<-tmp[,2]
    
    datPow$powerShift[j]<-sum(resInt)/B
    datPow$powerSlope[j]<-sum(resSlp)/B
    cat(paste(sep="","power to detect a shift in ab = ",datPow$powerShift[j],", power to detect a change in ab-age slope = ",datPow$powerSlope[j]," [",datPow$effectSize[j],", ",datPow$n[j],", ",datPow$pDropOut[j],"]\n\n"))
  }
  
  # save(list=ls(),file="Vacc-iNTS_SampleSize_MalariaInfection.RData")
  datPowMalnAb <- datPow
  export(datPow, here("data", "datPowMalnAb.RData"))
  
  # produce power curves
  pdf(file="~/work/PhD/VacciNTS/Sample size/Vacc-iNTS_SampleSize_MalariaInfection_powerCurves_pDO0.2.pdf",width=16,height=9)
  ggplot(data=datPowMalnAb[datPowMalnAb$pDropOut==0.2,],mapping=aes(x=effectSize,y=powerShift,group=as.factor(n),color=as.factor(n))) +
    geom_point(size=3) +
    geom_line(lwd=1.5) +
    geom_hline(yintercept=0.8,lwd=1,lty=2,col="darkgrey") + 
    scale_color_manual(values=colorRampPalette(brewer.pal(8,"Spectral"))(length(levels(factor(datPow$n)))),name="Number of paired samples") +
    scale_x_continuous(breaks = seq(0,50,1)) +
    ylab("Power") +
    xlab("Effect size (shift in average ab)") +
    ggtitle("Power curves: impact of Malnutrition on NTS antibody responses, assuming 20% drop-out.") +
    theme(text=element_text(size=16)) +
    theme_bw()
  dev.off()
  
  # produce a few example figures of data generated
  c1<-col2rgb("steelblue")
  c1<-rgb(c1[1],c1[2],c1[3],alpha=75,maxColorValue=255)
  c2<-col2rgb("orange")
  c2<-rgb(c2[1],c2[2],c2[3],alpha=75,maxColorValue=255)
  
  pdf(file="Vacc-iNTS_SampleSize_MalariaInfection_n1000_es0.1_pDO0.2.pdf",width=16,height=9)
  datSim<-simulateDataMalnAb(n=1000,prevBin=0.2,abSD=30,abSlopeMonth=2.5,iccSD=0.03,effectSize=-12.5,effectSizeSD=30,pDropOut=0.2,ageMin=8,ageMax=36)
  datSim$pid_n <- as.numeric(stringr::str_extract(datSim$pid, "\\d+"))
  mod<-gee(ab~agem8*exposure,id=pid_n,data=datSim,corstr="exchangeable")
  
  g<-ggplot(data=datSim,mapping=aes(x=agem8,y=ab,col=exposure, group = pid_n)) +
    geom_point() + 
    geom_line() +
    scale_color_manual(values=c(c1,c2)) +
    geom_abline(lwd=1,intercept=coef(mod)["(Intercept)"]+coef(mod)["exposureexposed"],slope=coef(mod)["agem8"]+coef(mod)["agem8:exposureexposed"],col="orange") + 
    geom_abline(lwd=1,intercept=coef(mod)["(Intercept)"],slope=coef(mod)["agem8"],col="steelblue") +
    theme_bw() +
    xlab("Age in months") +
    ylab("Antibody titre") +
    ggtitle("Example results for one simulation: antibody responses with age by malaria exposure")
  print(g)
  dev.off()
  
  
  ###############################################
  ## sickle cell
  
  if(is.na(args[2]) | args[2]=="sickle"){
    # simulate data
    simulateDataSSAb<-function(n,prevBin=0.01,abSD=30,abSlopeMonth=2.5,iccSD=0.03,effectSize=-12.5,effectSizeSD=30,pDropOut=0.2,ageMin=8,ageMax=36){
      # n<-100 # this would yield 2n measurements, each pair of measurements for the n individuals 3 months apart
      # prevBin<-0.2 # malaria prevalence in children 6-60 months, Southern Malawi, MIS 2017 p.43 Fig 4.4, https://dhsprogram.com/pubs/pdf/MIS28/MIS28.pdf
      # abSD<- 39.8 # Nyirenda et al (2014)
      # abSlopeMonth<-(3.2) # Nyirenda et al (2014)
      # ageMin<-8; ageMax<-28 # Nyirenda et al (2014)
      
      pid<-paste(sep="","p",1:n)
      ages<-runif(n=n,min=ageMin,max=ageMax) 
      binCov<-rbinom(n=2*n,size=1,prob=prevBin) # allows participants to change state for each visit; assumed independent risk at each visit and risk the same at all ages -- both of which may not be entirely realistic...
      
      u<-rnorm(n=n,mean=0,sd=iccSD)
      v<-rnorm(n=n,mean=0,sd=iccSD)
      ab8m<-rnorm(n=n,mean=0,sd=abSD)
      
      datSim<-data.frame(pid=rep(pid,2),
                         visit=c(rep(1,n),rep(2,n)),
                         age=c(ages,ages+rnorm(n,mean=3,sd=0.05)), # approximately 3 months between initial recruitment and follow-up but allowing minor deviations from planned visit dates
                         ageSlopeEffectLimit=NA,
                         exposure=ifelse(binCov==1,"exposed","unexposed"),
                         u=rep(u,2),
                         v=rep(v,2),
                         ab8m=rep(ab8m,2),
                         ab=NA)
      
      datSim$exposure<-factor(datSim$exposure,levels=c("unexposed","exposed"))
      datSim$agem8<-datSim$age-8
      
      for(i in 1:nrow(datSim)){
        pid<-datSim$pid[i]
        idx1<-which(datSim$pid==pid & datSim$visit==1)
        idx2<-which(datSim$pid==pid & datSim$visit==2)
        #datSim$ageSlopeEffectLimit[i]<-ifelse(datSim$exposure[idx1]=="exposed",datSim$age[idx1]-0.5,datSim$age[i])
        datSim$ageSlopeEffectLimit[idx2]<-ifelse(datSim$exposure[idx1]=="exposed",datSim$age[idx2]-0.5,datSim$age[idx2])
      }
      
      #datSim$ab<-datSim$ab8m+datSim$u+(abSlopeMonth+datSim$v)*(datSim$age-8-ifelse(datSim$visit==1 & datSim$exposure=="exposed",1,ifelse(datSim$visit==2 & datSim$exposure=="exposed",0.5,0)))+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
      datSim$ab<-datSim$ab8m+datSim$u+(abSlopeMonth+datSim$v)*(datSim$ageSlopeEffectLimit-8)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
      #datSim$ab<-datSim$ab8m+datSim$u+(abSlopeMonth+datSim$v)*(datSim$age-8)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
      datSim$ab[datSim$visit==2 & rbinom(n=n,size=1,prob=pDropOut)==1]<-NA
      datSim$ab[datSim$ab<=0]<-0.0001
      
      return(datSim)
    }
    
    # detect a change in slope of ab with age for i) malaria cases, ii) malaria vaccinated, iii) other binary covariates (linear mixed model approach)
    fitGEETestCoeff<-function(dat,alpha=0.05,coefName1="exposureexposed",coefName2=NULL){
      library(gee)
      #mod<-lmer(ab~agem8*exposure+(1|pid),REML=F,data=dat) # REML=F as we test fixed parameters; cannot have a second random effect even though this is how data was generated; would need more data points per participant...
      mod<-gee(ab~agem8*exposure,id=as.factor(pid),data=dat[order(dat$pid,dat$visit),],corstr="exchangeable",na.action=na.omit)
      p1<-2*(1-pnorm(abs(summary(mod)$coefficients[coefName1,"Robust z"])))
      p2<-2*(1-pnorm(abs(summary(mod)$coefficients[coefName2,"Robust z"])))
      
      #if(summary(mod)$coefficients[coefName1,"Pr(>|t|)"]<alpha){
      if(p1<alpha){
        res<-1
      }else{
        res<-0
      }
      
      if(!is.null(coefName2)){
        #if(summary(mod)$coefficients[coefName2,"Pr(>|t|)"]<alpha){
        if(p2<alpha){
          res<-c(res,1)
        }else{
          res<-c(res,0)
        }
      }
      
      return(res)
    }
  }
  
  
  # Sickle cell prevalence ~2%; try sample sizes 2000 and effect sizes 0.1, 0.2, 0.3
  # gr<-expand.grid(c(2500),c(12),c(0.2))
  gr<-expand.grid(c(2500),c(50,60,70,80),c(0.2))
  #gr<-expand.grid(c(2000, 2500, 3000),c(30,40,50),c(0.2))
  datPow<-data.frame(
    alpha=0.05,
    sicklePrev=0.01,
    effectSize=gr[,2],
    n=gr[,1],
    pDropOut=gr[,3],
    powerShift=NA,
    powerSlope=NA
  )
  
  for(j in 1:nrow(datPow)){
    cat(paste(sep="","## effectSize = ",datPow$effectSize[j],", n = ",datPow$n[j],", pDropOut = ",datPow$pDropOut[j]," ##\n"))
    
    #resInt<-rep(NA,B)
    #resSlp<-rep(NA,B)
    
    res<-foreach(i = 1:B) %dopar% {
      #if(i %% 100 == 0){
      #	cat(paste(sep="","..",i,"\n"))
      #}
      
      datSim<-simulateDataSSAb(n=datPow$n[j],prevBin=datPow$sicklePrev[j],abSD=30,abSlopeMonth=(2.5),iccSD=0.03,effectSize=datPow$effectSize[j],effectSizeSD=30,pDropOut=datPow$pDropOut[j],ageMin=8,ageMax=36)
      tmp<-fitGEETestCoeff(dat=datSim,coefName1="exposureexposed",coefName2="agem8:exposureexposed",alpha=datPow$alpha[j])
      #resInt[i]<-tmp[1]
      #resSlp[i]<-tmp[2]
      tmp
    }
    
    tmp<-matrix(unlist(res),byrow=T,ncol=2)
    resInt<-tmp[,1]
    resSlp<-tmp[,2]
    
    datPow$powerShift[j]<-sum(resInt)/B
    datPow$powerSlope[j]<-sum(resSlp)/B
    cat(paste(sep="","power to detect a shift in ab = ",datPow$powerShift[j],", power to detect a change in ab-age slope = ",datPow$powerSlope[j]," [",datPow$effectSize[j],", ",datPow$n[j],", ",datPow$pDropOut[j],"]\n\n"))
  }
  
  # save(list=ls(),file="Vacc-iNTS_SampleSize_MalariaInfection.RData")
  datPowSickleAb <- datPow
  export(datPowSickleAb, here("data", "datPowSickleAb.RData"))
  
  # produce power curves
  pdf(file="~/work/PhD/VacciNTS/Sample size/Vacc-iNTS_SampleSize_sicklecell_powerCurves_pDO0.2.pdf",width=16,height=9)
  ggplot(data=datPowSickleAb[datPowSickleAb$pDropOut==0.2,],mapping=aes(x=effectSize,y=powerShift,group=as.factor(n),color=as.factor(n))) +
    geom_point(size=3) +
    geom_line(lwd=1.5) +
    geom_hline(yintercept=0.8,lwd=1,lty=2,col="darkgrey") + 
    scale_color_manual(values=colorRampPalette(brewer.pal(8,"Spectral"))(length(levels(factor(datPow$n)))),name="Number of paired samples") +
    scale_x_continuous(breaks = seq(0,80,1)) +
    ylab("Power") +
    xlab("Effect size (shift in average ab)") +
    ggtitle("Power curves: impact of sickle cell status on NTS antibody responses, assuming 20% drop-out.") +
    theme(text=element_text(size=16)) +
    theme_bw()
  dev.off()
  
  # produce a few example figures of data generated
  c1<-col2rgb("steelblue")
  c1<-rgb(c1[1],c1[2],c1[3],alpha=75,maxColorValue=255)
  c2<-col2rgb("orange")
  c2<-rgb(c2[1],c2[2],c2[3],alpha=75,maxColorValue=255)
  
  pdf(file="Vacc-iNTS_SampleSize_MalariaInfection_n1000_es0.1_pDO0.2.pdf",width=16,height=9)
  datSim<-simulateDataSSAb(n=1000,prevBin=0.2,abSD=30,abSlopeMonth=2.5,iccSD=0.03,effectSize=-12.5,effectSizeSD=30,pDropOut=0.2,ageMin=8,ageMax=36)
  datSim$pid_n <- as.numeric(stringr::str_extract(datSim$pid, "\\d+"))
  mod<-gee(ab~agem8*exposure,id=pid_n,data=datSim,corstr="exchangeable")
  
  g<-ggplot(data=datSim,mapping=aes(x=agem8,y=ab,col=exposure, group = pid_n)) +
    geom_point() + 
    geom_line() +
    scale_color_manual(values=c(c1,c2)) +
    geom_abline(lwd=1,intercept=coef(mod)["(Intercept)"]+coef(mod)["exposureexposed"],slope=coef(mod)["agem8"]+coef(mod)["agem8:exposureexposed"],col="orange") + 
    geom_abline(lwd=1,intercept=coef(mod)["(Intercept)"],slope=coef(mod)["agem8"],col="steelblue") +
    theme_bw() +
    xlab("Age in months") +
    ylab("Antibody titre") +
    ggtitle("Example results for one simulation: antibody responses with age by malaria exposure")
  print(g)
  dev.off()
  
  
##############################
##  Enteric infection
  
  
  if(is.na(args[2]) | args[2]=="enteric"){
    # simulate data
    simulateDataEnAb<-function(n,prevBin=0.08,abSD=30,abSlopeMonth=2.5,iccSD=0.03,effectSize=-12.5,effectSizeSD=30,pDropOut=0.2,ageMin=8,ageMax=36){
      # n<-100 # this would yield 2n measurements, each pair of measurements for the n individuals 3 months apart
      # prevBin<-0.2 # malaria prevalence in children 6-60 months, Southern Malawi, MIS 2017 p.43 Fig 4.4, https://dhsprogram.com/pubs/pdf/MIS28/MIS28.pdf
      # abSD<- 39.8 # Nyirenda et al (2014)
      # abSlopeMonth<-(3.2) # Nyirenda et al (2014)
      # ageMin<-8; ageMax<-28 # Nyirenda et al (2014)
      
      pid<-paste(sep="","p",1:n)
      ages<-runif(n=n,min=ageMin,max=ageMax) 
      binCov<-rbinom(n=2*n,size=1,prob=prevBin) # allows participants to change state for each visit; assumed independent risk at each visit and risk the same at all ages -- both of which may not be entirely realistic...
      
      u<-rnorm(n=n,mean=0,sd=iccSD)
      v<-rnorm(n=n,mean=0,sd=iccSD)
      ab8m<-rnorm(n=n,mean=0,sd=abSD)
      
      datSim<-data.frame(pid=rep(pid,2),
                         visit=c(rep(1,n),rep(2,n)),
                         age=c(ages,ages+rnorm(n,mean=3,sd=0.05)), # approximately 3 months between initial recruitment and follow-up but allowing minor deviations from planned visit dates
                         ageSlopeEffectLimit=NA,
                         exposure=ifelse(binCov==1,"exposed","unexposed"),
                         u=rep(u,2),
                         v=rep(v,2),
                         ab8m=rep(ab8m,2),
                         ab=NA)
      
      datSim$exposure<-factor(datSim$exposure,levels=c("unexposed","exposed"))
      datSim$agem8<-datSim$age-8
      
      for(i in 1:nrow(datSim)){
        pid<-datSim$pid[i]
        idx1<-which(datSim$pid==pid & datSim$visit==1)
        idx2<-which(datSim$pid==pid & datSim$visit==2)
        #datSim$ageSlopeEffectLimit[i]<-ifelse(datSim$exposure[idx1]=="exposed",datSim$age[idx1]-0.5,datSim$age[i])
        datSim$ageSlopeEffectLimit[idx2]<-ifelse(datSim$exposure[idx1]=="exposed",datSim$age[idx2]-0.5,datSim$age[idx2])
       # datSim$ageSlopeEffectLimit[idx1]<-datSim$age[idx1]
        
      }
      
      #datSim$ab<-datSim$ab8m+datSim$u+(abSlopeMonth+datSim$v)*(datSim$age-8-ifelse(datSim$visit==1 & datSim$exposure=="exposed",1,ifelse(datSim$visit==2 & datSim$exposure=="exposed",0.5,0)))+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
      datSim$ab<-datSim$ab8m+datSim$u+(abSlopeMonth+datSim$v)*(datSim$ageSlopeEffectLimit-8)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
      #datSim$ab<-datSim$ab8m+datSim$u+(abSlopeMonth+datSim$v)*(datSim$age-8)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
      datSim$ab[datSim$visit==2 & rbinom(n=n,size=1,prob=pDropOut)==1]<-NA
      datSim$ab[datSim$ab<=0]<-0.0001
      
      return(datSim)
    }
    
    # detect a change in slope of ab with age for i) malaria cases, ii) malaria vaccinated, iii) other binary covariates (linear mixed model approach)
    fitGEETestCoeff<-function(dat,alpha=0.05,coefName1="exposureexposed",coefName2=NULL){
      library(gee)
      #mod<-lmer(ab~agem8*exposure+(1|pid),REML=F,data=dat) # REML=F as we test fixed parameters; cannot have a second random effect even though this is how data was generated; would need more data points per participant...
      mod<-gee(ab~agem8*exposure,id=as.factor(pid),data=dat[order(dat$pid,dat$visit),],corstr="exchangeable",na.action=na.omit)
      p1<-2*(1-pnorm(abs(summary(mod)$coefficients[coefName1,"Robust z"])))
      p2<-2*(1-pnorm(abs(summary(mod)$coefficients[coefName2,"Robust z"])))
      
      #if(summary(mod)$coefficients[coefName1,"Pr(>|t|)"]<alpha){
      if(p1<alpha){
        res<-1
      }else{
        res<-0
      }
      
      if(!is.null(coefName2)){
        #if(summary(mod)$coefficients[coefName2,"Pr(>|t|)"]<alpha){
        if(p2<alpha){
          res<-c(res,1)
        }else{
          res<-c(res,0)
        }
      }
      
      return(res)
    }
  }
  
  
  # enteric infection; prevalence ~8%; try sample sizes 2000 and effect sizes 0.1, 0.2, 0.3
  gr<-expand.grid(c(2000, 2500, 3000),c(12,15,18,22,25,30),c(0.2))
  #gr<-expand.grid(c(2500),c(12),c(0.2))
  datPow<-data.frame(
    alpha=0.05,
    malariaPrev=0.16,
    effectSize=gr[,2],
    n=gr[,1],
    pDropOut=gr[,3],
    powerShift=NA,
    powerSlope=NA
  )
  
  for(j in 1:nrow(datPow)){
    cat(paste(sep="","## effectSize = ",datPow$effectSize[j],", n = ",datPow$n[j],", pDropOut = ",datPow$pDropOut[j]," ##\n"))
    
    #resInt<-rep(NA,B)
    #resSlp<-rep(NA,B)
    
    res<-foreach(i = 1:B) %dopar% {
      #if(i %% 100 == 0){
      #	cat(paste(sep="","..",i,"\n"))
      #}
      
      datSim<-simulateDataEnAb(n=datPow$n[j],prevBin=datPow$malariaPrev[j],abSD=30,abSlopeMonth=(2.5),iccSD=0.03,effectSize=datPow$effectSize[j],effectSizeSD=30,pDropOut=datPow$pDropOut[j],ageMin=8,ageMax=36)
      tmp<-fitGEETestCoeff(dat=datSim,coefName1="exposureexposed",coefName2="agem8:exposureexposed",alpha=datPow$alpha[j])
      #resInt[i]<-tmp[1]
      #resSlp[i]<-tmp[2]
      tmp
    }
    
    tmp<-matrix(unlist(res),byrow=T,ncol=2)
    resInt<-tmp[,1]
    resSlp<-tmp[,2]
    
    datPow$powerShift[j]<-sum(resInt)/B
    datPow$powerSlope[j]<-sum(resSlp)/B
    cat(paste(sep="","power to detect a shift in ab = ",datPow$powerShift[j],", power to detect a change in ab-age slope = ",datPow$powerSlope[j]," [",datPow$effectSize[j],", ",datPow$n[j],", ",datPow$pDropOut[j],"]\n\n"))
  }
  
  datPowEntericInfAb <- datPow
  export(datPow, here("data", "datPowEntericInfAb.RData"))
  
  # produce power curves
  pdf(file="~/work/PhD/VacciNTS/Sample size/Vacc-iNTS_SampleSize_datPowEntericInfAb.pdf",width=16,height=9)
  ggplot(data=datPowEntericInfAb[datPowEntericInfAb$pDropOut==0.2,],mapping=aes(x=effectSize,y=powerShift,group=as.factor(n),color=as.factor(n))) +
    geom_point(size=3) +
    geom_line(lwd=1.5) +
    geom_hline(yintercept=0.8,lwd=1,lty=2,col="darkgrey") + 
    scale_color_manual(values=colorRampPalette(brewer.pal(8,"Spectral"))(length(levels(factor(datPow$n)))),name="Number of paired samples") +
    scale_x_continuous(breaks = seq(0,20,1)) +
    ylab("Power") +
    xlab("Effect size (shift in average ab)") +
    ggtitle("Power curves: Impact of enteric exposure on NTS antibody responses, drop-out (20%)") +
    theme(text=element_text(size=16)) +
    theme_bw()
  dev.off()