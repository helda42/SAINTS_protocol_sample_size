# Sample size calculation for Vacc-iNTS

rm(list = ls())

library(ggplot2)
library(RColorBrewer)
#library(lme4)
#library(lmerTest)
library(gee) # using GEE instead of LMM; fixed intercept in LMM captured some of the effect due to exposure...
library(doParallel)
library(pwr)


########################
## 1. Malaria exposure
#######################

args<-commandArgs(T)

if(!is.na(args[1])){
  cat(paste(sep="","Using ",args[1]," cores!\n\n"))
  #cl <- makeCluster(args[1])
  #registerDoParallel(cl)
  registerDoParallel(cores=args[1])
}else{
  cat(paste(sep="","Using ",2," cores!\n\n"))
  #cl <- makeCluster(2)
  #registerDoParallel(cl)
  registerDoParallel(cores=2)
}

set.seed(20200128)

# power parameters
alpha<-0.05
power<-0.8
B<-1e3 # number of simulations

#setwd("~/work/MLW_LSTM/DaleHelen_StatsSupport/2021-10-20_SimulationScript/output")


if(is.na(args[2]) | args[2]=="malaria"){
  # simulate data
  simulateData<-function(n,prevBin=0.14,sbaSD=0.4,sbaSlopeMonth=(-0.1),iccSD=0.03,effectSize=0.5,effectSizeSD=0.25,pDropOut=0.2,ageMin=9,ageMax=36){
    # n<-100 # this would yield 2n measurements, each pair of measurements for the n individuals 3 months apart
    # prevBin<-0.26 # malaria prevalence in children 6-60 months, Southern Malawi, MIS 2017 p.43 Fig 4.4, https://dhsprogram.com/pubs/pdf/MIS28/MIS28.pdf
    # sbaSD<-0.4 # Nyirenda et al (2014)
    # sbaSlopeMonth<-(-0.105) # Nyirenda et al (2014)
    # ageMin<-9; ageMax<-36 # Nyirenda et al (2014)
    
    pid<-paste(sep="","p",1:n)
    ages<-runif(n=n,min=ageMin,max=ageMax) 
    binCov<-rbinom(n=2*n,size=1,prob=prevBin) # allows participants to change state for each visit; assumed independent risk at each visit and risk the same at all ages -- both of which may not be entirely realistic...
    
    u<-rnorm(n=n,mean=0,sd=iccSD)
    v<-rnorm(n=n,mean=0,sd=iccSD)
    sba9m<-rnorm(n=n,mean=0,sd=sbaSD)
    
    datSim<-data.frame(pid=rep(pid,2),
                       visit=c(rep(1,n),rep(2,n)),
                       age=c(ages,ages+rnorm(n,mean=3,sd=0.05)), # approximately 3 months between initial recruitment and follow-up but allowing minor deviations from planned visit dates
                       ageSlopeEffectLimit=NA,
                       exposure=ifelse(binCov==1,"exposed","unexposed"),
                       u=rep(u,2),
                       v=rep(v,2),
                       sba9m=rep(sba9m,2),
                       sba=NA)
    
    datSim$exposure<-factor(datSim$exposure,levels=c("unexposed","exposed"))
    datSim$agem9<-datSim$age-9
    
    datSim$ageSlopeEffectLimit[i]<-datSim$age[i]
    for(i in 1:nrow(datSim)){
      pid<-datSim$pid[i]
      idx1<-which(datSim$pid==pid & datSim$visit==1)
      idx2<-which(datSim$pid==pid & datSim$visit==2)
      datSim$ageSlopeEffectLimit[idx2]<-ifelse(datSim$exposure[idx1]=="exposed",datSim$age[idx2]-0.5,datSim$age[idx2])
    }
    
    #datSim$sba<-datSim$sba9m+datSim$u+(sbaSlopeMonth+datSim$v)*(datSim$age-9-ifelse(datSim$visit==1 & datSim$exposure=="exposed",1,ifelse(datSim$visit==2 & datSim$exposure=="exposed",0.5,0)))+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
    datSim$sba<-datSim$sba9m+datSim$u+(sbaSlopeMonth+datSim$v)*(datSim$ageSlopeEffectLimit-9)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
    #datSim$sba<-datSim$sba9m+datSim$u+(sbaSlopeMonth+datSim$v)*(datSim$age-9)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
    datSim$sba[datSim$visit==2 & rbinom(n=n,size=1,prob=pDropOut)==1]<-NA
    
    return(datSim)
  }
  
  # detect a change in slope of SBA with age for i) malaria cases, ii) malaria vaccinated, iii) other binary covariates (linear mixed model approach)
  fitGEETestCoeff<-function(dat,alpha=0.05,coefName1="exposureexposed",coefName2=NULL){
    library(gee)
    #mod<-lmer(sba~agem9*exposure+(1|pid),REML=F,data=dat) # REML=F as we test fixed parameters; cannot have a second random effect even though this is how data was generated; would need more data points per participant...
    mod<-gee(sba~agem9*exposure,id=as.factor(pid),data=dat[order(dat$pid,dat$visit),],corstr="exchangeable",na.action=na.omit)
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
  
  
  # malaria infection; prevalence ~14%; try sample sizes 1000, 1500, 2000 and effect sizes 0.05, 0.1, 0.15
  gr<-expand.grid(c(2000, 2500, 3000),c(0.2, 0.225, 0.25, 0.275, 0.3),c(0.2))
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
      
      datSim<-simulateData(n=datPow$n[j],prevBin=datPow$malariaPrev[j],sbaSD=0.4,sbaSlopeMonth=(-0.1),iccSD=0.03,effectSize=datPow$effectSize[j],effectSizeSD=0.25,pDropOut=datPow$pDropOut[j],ageMin=9,ageMax=36)
      tmp<-fitGEETestCoeff(dat=datSim,coefName1="exposureexposed",coefName2="agem9:exposureexposed",alpha=datPow$alpha[j])
      #resInt[i]<-tmp[1]
      #resSlp[i]<-tmp[2]
      tmp
    }
    
    tmp<-matrix(unlist(res),byrow=T,ncol=2)
    resInt<-tmp[,1]
    resSlp<-tmp[,2]
    
    datPow$powerShift[j]<-sum(resInt)/B
    datPow$powerSlope[j]<-sum(resSlp)/B
    cat(paste(sep="","power to detect a shift in SBA = ",datPow$powerShift[j],", power to detect a change in SBA-age slope = ",datPow$powerSlope[j]," [",datPow$effectSize[j],", ",datPow$n[j],", ",datPow$pDropOut[j],"]\n\n"))
  }
}

#stopCluster(cl)

save(list=ls(),file="~/work/PhD/VacciNTS/Sample size/SAINTS_protocol_sample_size/data/Vacc-iNTS_SampleSize_MalariaInfection.RData")


# produce power curves
pdf(file="~/work/PhD/VacciNTS/Sample size/SAINTS_protocol_sample_size/output/Vacc-iNTS_SampleSize_MalariaInfection_powerCurves_pDO0.2.pdf",width=16,height=9)
ggplot(data=datPow[datPow$pDropOut==0.2,],mapping=aes(x=effectSize,y=powerShift,group=as.factor(n),color=as.factor(n))) +
  geom_point(size=3) +
  geom_line(lwd=1.5) +
  geom_hline(yintercept=0.8,lwd=2,lty=2,col="darkgrey") + 
  scale_color_manual(values=colorRampPalette(brewer.pal(9,"Spectral"))(length(levels(factor(datPow$n)))),name="number of paired samples") +
  ylab("power") +
  xlab("effect size (shift in average SBA)") +
  ggtitle("Power for different effect sizes and sample sizes, assuming 20% drop-out.") +
  theme(text=element_text(size=16))
dev.off()

# produce a few example figures of data generated
c1<-col2rgb("steelblue")
c1<-rgb(c1[1],c1[2],c1[3],alpha=75,maxColorValue=255)
c2<-col2rgb("orange")
c2<-rgb(c2[1],c2[2],c2[3],alpha=75,maxColorValue=255)

pdf(file="~/work/PhD/VacciNTS/Sample size/SAINTS_protocol_sample_size/output/Vacc-iNTS_SampleSize_MalariaInfection_n1000_es0.1_pDO0.2.pdf",width=16,height=9)
datSim<-simulateData(n=1000,prevBin=0.14,sbaSD=0.4,sbaSlopeMonth=(-0.1),iccSD=0.03,effectSize=0.1,effectSizeSD=0.25,pDropOut=0.2,ageMin=9,ageMax=36)
mod<-gee(sba~agem9*exposure,id=as.factor(pid),data=datSim[order(datSim$pid,datSim$visit),],corstr="exchangeable",na.action=na.omit)

g<-ggplot(data=datSim,mapping=aes(x=agem9,y=sba,col=exposure,group=pid)) +
  geom_point() + 
  geom_line() +
  scale_color_manual(values=c(c1,c2)) +
  geom_abline(lwd=2,intercept=coef(mod)["(Intercept)"]+coef(mod)["exposureexposed"],slope=coef(mod)["agem9"]+coef(mod)["agem9:exposureexposed"],col="orange") + 
  geom_abline(lwd=2,intercept=coef(mod)["(Intercept)"],slope=coef(mod)["agem9"],col="steelblue")
print(g)
dev.off()

pdf(file="~/work/PhD/VacciNTS/Sample size/SAINTS_protocol_sample_size/output/Vacc-iNTS_SampleSize_MalariaInfection_n2000_es0.1_pDO0.2.pdf",width=16,height=9)
datSim<-simulateData(n=2000,prevBin=0.14,sbaSD=0.4,sbaSlopeMonth=(-0.1),iccSD=0.03,effectSize=0.1,effectSizeSD=0.25,pDropOut=0.2,ageMin=9,ageMax=36)

g<-ggplot(data=datSim,mapping=aes(x=agem9,y=sba,col=exposure,group=pid)) +
  geom_point() + 
  geom_line() +
  scale_color_manual(values=c(c1,c2)) +
  geom_abline(lwd=2,slope=-0.1,intercept=0.1,col="orange") + 
  geom_abline(lwd=2,slope=-0.1,intercept=0,col="steelblue")
print(g)
dev.off()

pdf(file="~/work/PhD/VacciNTS/Sample size/SAINTS_protocol_sample_size/output/Vacc-iNTS_SampleSize_MalariaInfection_n2000_es0.15_pDO0.2.pdf",width=16,height=9)
datSim<-simulateData(n=2000,prevBin=0.14,sbaSD=0.4,sbaSlopeMonth=(-0.1),iccSD=0.03,effectSize=0.15,effectSizeSD=0.25,pDropOut=0.2,ageMin=9,ageMax=36)
mod<-gee(sba~agem9*exposure,id=as.factor(pid),data=datSim[order(datSim$pid,datSim$visit),],corstr="exchangeable",na.action=na.omit)

g<-ggplot(data=datSim,mapping=aes(x=agem9,y=sba,col=exposure,group=pid)) +
  geom_point() + 
  geom_line() +
  scale_color_manual(values=c(c1,c2)) +
  geom_abline(lwd=2,intercept=coef(mod)["(Intercept)"]+coef(mod)["exposureexposed"],slope=coef(mod)["agem9"]+coef(mod)["agem9:exposureexposed"],col="orange") + 
  geom_abline(lwd=2,intercept=coef(mod)["(Intercept)"],slope=coef(mod)["agem9"],col="steelblue")
print(g)
dev.off()

pdf(file="~/work/PhD/VacciNTS/Sample size/SAINTS_protocol_sample_size/output/Vacc-iNTS_SampleSize_MalariaInfection_n2000_es0.2_pDO0.2.pdf",width=16,height=9)
datSim<-simulateData(n=2000,prevBin=0.14,sbaSD=0.4,sbaSlopeMonth=(-0.1),iccSD=0.03,effectSize=0.2,effectSizeSD=0.25,pDropOut=0.2,ageMin=9,ageMax=36)
mod<-gee(sba~agem9*exposure,id=as.factor(pid),data=datSim[order(datSim$pid,datSim$visit),],corstr="exchangeable",na.action=na.omit)

g<-ggplot(data=datSim,mapping=aes(x=agem9,y=sba,col=exposure,group=pid)) +
  geom_point() + 
  geom_line() +
  scale_color_manual(values=c(c1,c2)) +
  geom_abline(lwd=2,intercept=coef(mod)["(Intercept)"]+coef(mod)["exposureexposed"],slope=coef(mod)["agem9"]+coef(mod)["agem9:exposureexposed"],col="orange") + 
  geom_abline(lwd=2,intercept=coef(mod)["(Intercept)"],slope=coef(mod)["agem9"],col="steelblue")
print(g)
dev.off()

##############################
## 2. Malaria vaccination ##
############################

if(is.na(args[2]) | args[2]=="malvacc"){
  
  simulateDataMalVacc<-function(n,prevBin=0.26,dirVE=0.4,indirVE=0.1,vaccCov=0.6,sbaSD=0.4,sbaSlopeMonth=(-0.1),iccSD=0.03,effectSize=0.5,effectSizeSD=0.25,pDropOut=0.1,ageMin=9,ageMax=36){
    # n<-100 # this would yield 2n measurements, each pair of measurements for the n individuals 3 months apart
    # prevBin<-0.26 # malaria prevalence in children 6-60 months, Southern Malawi, MIS 2017 p.43 Fig 4.4, https://dhsprogram.com/pubs/pdf/MIS28/MIS28.pdf
    # sbaSD<-0.4 # Nyirenda et al (2014)
    # sbaSlopeMonth<-(-0.105) # Nyirenda et al (2014)
    # ageMin<-9; ageMax<-36 # Nyirenda et al (2014)
    
    pid<-paste(sep="","p",1:n)
    ages<-runif(n=n,min=ageMin,max=ageMax) 
    binCov<-rbinom(n=2*n,size=1,prob=prevBin) # allows participants to change state for each visit; assumed independent risk at each visit and risk the same at all ages -- both of which may not be entirely realistic...
    
    u<-rnorm(n=n,mean=0,sd=iccSD)
    v<-rnorm(n=n,mean=0,sd=iccSD)
    sba9m<-rnorm(n=n,mean=0,sd=sbaSD)
    
    cluster<-c(rep(c(rep(0,n/2),rep(1,n/2)),2)) # vaccine cluster membership
    vaccTmp<-sample(0:1,size=n/2,replace=T,prob=c(1-vaccCov,vaccCov))
    vaccine<-c(rep(0,n/2),vaccTmp,rep(0,n/2),vaccTmp) # vaccination recipients
    
    protCov<-rep(0,2*n)
    idxIndi<-which(binCov==1 & cluster==1 & vaccine==0)
    idxVacc<-which(binCov==1 & vaccine==1)
    protCov[sample(idxIndi,replace=T,size=round(indirVE*length(idxIndi)))]<-1
    protCov[sample(idxVacc,replace=T,size=round((1-(1-indirVE)*(1-dirVE))*length(idxVacc)))]<-1
    
    datSim<-data.frame(pid=rep(pid,2),
                       visit=c(rep(1,n),rep(2,n)),
                       age=c(ages,ages+rnorm(n,mean=3,sd=0.05)), # approximately 3 months between initial recruitment and follow-up but allowing minor deviations from planned visit dates
                       vaccCluster=cluster,
                       vaccine=vaccine,
                       ageSlopeEffectLimit=NA,
                       exposure=ifelse(binCov==1,"exposed","unexposed"),
                       infection=ifelse(binCov==1 & protCov==0,"infected","uninfected"),
                       u=rep(u,2),
                       v=rep(v,2),
                       sba9m=rep(sba9m,2),
                       sba=NA)
    
    datSim$exposure<-factor(datSim$exposure,levels=c("unexposed","exposed"))
    datSim$infection<-factor(datSim$infection,levels=c("uninfected","infected"))
    datSim$agem9<-datSim$age-9
    
    for(i in 1:nrow(datSim)){
      pid<-datSim$pid[i]
      idx1<-which(datSim$pid==pid & datSim$visit==1)
      idx2<-which(datSim$pid==pid & datSim$visit==2)
      datSim$ageSlopeEffectLimit[i]<-ifelse(datSim$infection[idx1]=="infected",datSim$age[idx1]-0.5,datSim$age[i])
    }
    
    #datSim$sba<-datSim$sba9m+datSim$u+(sbaSlopeMonth+datSim$v)*(datSim$age-9-ifelse(datSim$visit==1 & datSim$exposure=="exposed",1,ifelse(datSim$visit==2 & datSim$exposure=="exposed",0.5,0)))+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
    datSim$sba<-datSim$sba9m+datSim$u+(sbaSlopeMonth+datSim$v)*(datSim$ageSlopeEffectLimit-9)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$infection=="infected",1,0)
    #datSim$sba<-datSim$sba9m+datSim$u+(sbaSlopeMonth+datSim$v)*(datSim$age-9)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
    datSim$sba[datSim$visit==2 & rbinom(n=n,size=1,prob=pDropOut)==1]<-NA
    
    return(datSim)
  }
  
  fitGEETestCoeffMalVacc<-function(dat,alpha=0.05,coefName1="vaccine",coefName2=NULL){
    library(gee)
    #mod<-lmer(sba~agem9*exposure+(1|pid),REML=F,data=dat) # REML=F as we test fixed parameters; cannot have a second random effect even though this is how data was generated; would need more data points per participant...
    mod<-gee(sba~agem9*exposure+vaccCluster+vaccine,id=pid,data=dat,corstr="exchangeable")
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
  
  gr<-expand.grid(c(2500),c(0.42,0.46,0.5),c(0.2),c(0.25,0.5,0.75,1))
  datPowMalVacc<-data.frame(
    alpha=0.05,
    malariaPrev=0.16,
    effectSize=gr[,4],
    dirVE=gr[,2],
    indirVE=0.1,
    vaccCov=0.6,
    n=gr[,1],
    pDropOut=gr[,3],
    powerDir=NA,
    powerIndir=NA
  )
  
  for(j in 1:nrow(datPowMalVacc)){
    cat(paste(sep="","## direct VE = ",datPowMalVacc$dirVE[j],", n = ",datPowMalVacc$n[j],", pDropOut = ",datPowMalVacc$pDropOut[j],", effectSize = ",datPowMalVacc$effectSize[j]," ##\n"))
    
    # resInt<-rep(NA,B)
    # resSlp<-rep(NA,B)
    
    res<-foreach(i = 1:B) %dopar% {
      #if(i %% 100 == 0){
      #	cat(paste(sep="","..",i,"\n"))
      #}
      
      datSim<-simulateDataMalVacc(n=datPowMalVacc$n[j],prevBin=datPowMalVacc$malariaPrev[j],sbaSD=0.4,sbaSlopeMonth=(-0.1),iccSD=0.03,effectSize=datPowMalVacc$effectSize[j],effectSizeSD=0.25,pDropOut=datPowMalVacc$pDropOut[j],ageMin=9,ageMax=36,dirVE=datPowMalVacc$dirVE[j],indirVE=datPowMalVacc$indirVE[j],vaccCov=datPowMalVacc$vaccCov[j])
      tmp<-fitGEETestCoeffMalVacc(dat=datSim,coefName1="vaccine",coefName2="vaccCluster",alpha=datPowMalVacc$alpha[j])
      # resInt[i]<-tmp[1]
      # resSlp[i]<-tmp[2]
      tmp
    }
    
    tmp<-matrix(unlist(res),byrow=T,ncol=2)
    resInt<-tmp[,1]
    resSlp<-tmp[,2]
    
    datPowMalVacc$powerDir[j]<-sum(resInt)/B
    datPowMalVacc$powerIndir[j]<-sum(resSlp)/B
    cat(paste(sep="","power to detect vaccine direct effect = ",datPowMalVacc$powerDir[j],", power to detect vaccine indirect effect = ",datPowMalVacc$powerIndir[j]," [",datPowMalVacc$effectSize[j],", ",datPowMalVacc$dirVE[j],", ",datPowMalVacc$n[j],", ",datPowMalVacc$pDropOut[j],"]\n\n"))
  }
  
  save(list=ls(),file="Vacc-iNTS_SampleSize_MalariaInfection.RData")
  
}

#dat<-read.csv("vacciNTS_Power_MalVacc.csv")
datPowMalVacc<-read.csv("Vacc-iNTS_SampleSize_MalariaInfection_20200311_fullTable_1000.csv")

g<-ggplot(data=datPowMalVacc,mapping=aes(x=effectSize,y=power,group=factor(dirVE),color=factor(dirVE))) +
  geom_point(size=3) +
  geom_line(lwd=1.25) +
  geom_hline(yintercept=0.8,lwd=2,lty=2,col="darkgrey") + 
  scale_color_manual(values=c("steelblue","orange"),name="Malaria vaccine efficacy.") +
  ylab("power") +
  xlab("effect size of malaria infection on SBA") +
  ggtitle("Power to detect an association between malaria vaccination and iNTS serology (10% drop-out rate).") +
  theme(text=element_text(size=16)) + 
  ylim(c(0,1)) +
  theme_bw()

print(g)

################################
## 3. Salmonella stool exposure
################################

gr<-expand.grid(seq(1000,2500,by=100),c(0.1,0.15,0.2,0.25),c(0.2,0.3,0.4,0.5),c(1/24,0.06,0.08,0.099))

datPowSeroConv<-data.frame(
  alpha=0.05,
  prev=gr[,4],
  pDropOut=0.2,
  nTot=gr[,1],
  nTotComp=NA,
  n1=NA,
  n2=NA,
  p1=gr[,2],
  p2=gr[,3],
  power=NA
)

datPowSeroConv$nTotComp<-floor((1-datPowSeroConv$pDropOut)*datPowSeroConv$nTot)
datPowSeroConv<-datPowSeroConv[datPowSeroConv$p2/datPowSeroConv$p1<3 & datPowSeroConv$p2/datPowSeroConv$p1>1,]

datPowSeroConv$n1<-round(datPowSeroConv$prev*datPowSeroConv$nTotComp)
datPowSeroConv$n2<-round((1-datPowSeroConv$prev)*datPowSeroConv$nTotComp)
datPowSeroConv$p1p2<-factor(paste(sep="","p1 = ",datPowSeroConv$p1,", p2 = ",datPowSeroConv$p2))

for(j in 1:nrow(datPowSeroConv)){
  datPowSeroConv$power[j]<-pwr.2p2n.test(h=ES.h(p1=datPowSeroConv$p1[j],p2=datPowSeroConv$p2[j]),n1=datPowSeroConv$n1[j],n2=datPowSeroConv$n2[j])$power
}

pdf(file="~/work/PhD/VacciNTS/Sample size/SAINTS_protocol_sample_size/outputs/Vacc-iNTS_SampleSize_SalmonellaStoolExposure_SBA8.pdf",width=16,height=9)
ggplot(data=datPowSeroConv[datPowSeroConv$prev==0.06,],mapping=aes(x=nTot,y=power,group=p1p2,color=p1p2)) +
  geom_point(size=4) +
  geom_line(lwd=2) +
  scale_color_manual(values=colorRampPalette(brewer.pal(9,"Spectral"))(length(levels(datPowSeroConv$p1p2))),name="prop. seroconverted in each group") +
  geom_hline(yintercept=0.8,lwd=2,lty=2,col="darkgrey") + 
  geom_vline(xintercept=2000,lwd=2,lty=2,col="darkgrey") +
  theme(legend.position = "right",text=element_text(size=16)) +
  xlab("number of paired samples") +
  ylab("power") +
  theme_bw() +
  ggtitle("Power curves for different combinations of 3-month seroconversion rates. Prevalence of exposure = 6%, drop-out = 20%.")
dev.off()
