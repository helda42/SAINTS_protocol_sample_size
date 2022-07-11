##############################################################
## Sample size for seroconversion between risk factor groups ##
##############################################################

rm(list = ls())

library(pwr)
library(ggplot2)
library(rio)
library(here)
#tinytex::install_tinytex()


## Malnutrition

# grid of sample size, effect size (proportion seroconverted in each group), standard deviation, prevalance of exposure
#  gr<-expand.grid(seq(1000,2500,by=100),c(0.1,0.2),c(0.4,0.6),c(1/24,0.074,0.08,0.099))
gr<-expand.grid(c(2000, 2225, 2500, 2750, 3000),c(0.15,0.2),c(0.4,0.6),c(0.02))

datPowmaln_sero<-data.frame(
  alpha=0.05,
  prev=gr[,4],
  pDropOut=0.2,
  nTot=gr[,1],
  nTotComp=NA,
  n1=NA,
  n2=NA,
  delta=gr[,2],
  sd=gr[,3],
  es=NA,
  power=NA
)

datPowmaln_sero$nTotComp<-floor((1-datPowmaln_sero$pDropOut)*datPowmaln_sero$nTot)
datPowmaln_sero$n1<-round(datPowmaln_sero$prev*datPowmaln_sero$nTotComp)
datPowmaln_sero$n2<-round((1-datPowmaln_sero$prev)*datPowmaln_sero$nTotComp)
datPowmaln_sero$es<-round(digits=2,datPowmaln_sero$delta/datPowmaln_sero$sd)

for(j in 1:nrow(datPowmaln_sero)){
  datPowmaln_sero$power[j]<-pwr.t2n.test(n1=datPowmaln_sero$n1[j],n2=datPowmaln_sero$n2[j],d=datPowmaln_sero$es[j],sig.level=datPowmaln_sero$alpha[j])$power
}

datPowmaln_sero$EsSd<-factor(paste(sep="","ES = ",datPowmaln_sero$delta,", SD = ",datPowmaln_sero$sd))

datPowmaln_sero$EsSd


#pdf(file="~/work/PhD/VacciNTS/Sample size/SAINTS_protocol_sample_size/outputs/Vacc-iNTS_SampleSize_malnutrition_seroconvert.pdf",width=16,height=9)
ggplot(data=datPowmaln_sero,mapping=aes(x=nTot,y=power,group=EsSd, colour=EsSd)) +
  geom_point(size=3) +
  geom_line(lwd=1.5) +
  scale_color_manual(values=c("steelblue","greenyellow","salmon","orange"),name="effect size and std. dev.") +
  geom_hline(yintercept=0.8,lwd=1,lty=2,col="darkgrey") + 
  theme(legend.position = "right",text=element_text(size=16)) +
  xlab("Number of paired samples") +
  ylab("Power") +
  ggtitle("Power curves; Seroconversion by malnutrition status (prev 2%), drop-out (20%)") +
  theme_bw()
#dev.off()


## Malaria

# grid of sample size, effect size (proportion seroconverted in each group), standard deviation, prevalance of exposure
#  gr<-expand.grid(seq(1000,2500,by=100),c(0.1,0.2),c(0.4,0.6),c(1/24,0.074,0.08,0.099))
gr<-expand.grid(c(2000, 2225, 2500, 2750, 3000),c(0.05, 0.1),c(0.4,0.6),c(0.16))

datPowMalInf_sero<-data.frame(
  alpha=0.05,
  prev=gr[,4],
  pDropOut=0.2,
  nTot=gr[,1],
  nTotComp=NA,
  n1=NA,
  n2=NA,
  delta=gr[,2],
  sd=gr[,3],
  es=NA,
  power=NA
)

datPowMalInf_sero$nTotComp<-floor((1-datPowMalInf_sero$pDropOut)*datPowMalInf_sero$nTot)
datPowMalInf_sero$n1<-round(datPowMalInf_sero$prev*datPowMalInf_sero$nTotComp)
datPowMalInf_sero$n2<-round((1-datPowMalInf_sero$prev)*datPowMalInf_sero$nTotComp)
datPowMalInf_sero$es<-round(digits=2,datPowMalInf_sero$delta/datPowMalInf_sero$sd)

for(j in 1:nrow(datPowMalInf_sero)){
  datPowMalInf_sero$power[j]<-pwr.t2n.test(n1=datPowMalInf_sero$n1[j],n2=datPowMalInf_sero$n2[j],d=datPowMalInf_sero$es[j],sig.level=datPowMalInf_sero$alpha[j])$power
}

datPowMalInf_sero$EsSd<-factor(paste(sep="","ES = ",datPowMalInf_sero$delta,", SD = ",datPowMalInf_sero$sd))

datPowMalInf_sero$EsSd
datPowMalInf_sero


pdf(file="~/work/PhD/VacciNTS/Sample size/SAINTS_protocol_sample_size/outputs/Vacc-iNTS_SampleSize_malaria_seroconvert.pdf",width=16,height=9)
ggplot(data=datPowMalInf_sero,mapping=aes(x=nTot,y=power,group=EsSd, colour=EsSd)) +
  geom_point(size=3) +
  geom_line(lwd=1.5) +
  scale_color_manual(values=c("steelblue","greenyellow","salmon","orange"),name="effect size and std. dev.") +
  geom_hline(yintercept=0.8,lwd=1,lty=2,col="darkgrey") + 
  theme(legend.position = "right",text=element_text(size=16)) +
  xlab("Number of paired samples") +
  ylab("Power") +
  theme_bw() +
  ggtitle("Power curves; Seroconversion by malaria status (prev 16%), drop-out (20%)") +
  theme_bw()
dev.off()


## Anaemia
# grid of sample size, effect size (proportion seroconverted in each group), standard deviation, prevalance of exposure
#  gr<-expand.grid(seq(1000,2500,by=100),c(0.1,0.2),c(0.4,0.6),c(1/24,0.074,0.08,0.099))
gr<-expand.grid(c(2000, 2225, 2500, 2750, 3000),c(0.05, 0.1),c(0.4,0.6),c(0.34))

datPowAnaem_sero<-data.frame(
  alpha=0.05,
  prev=gr[,4],
  pDropOut=0.2,
  nTot=gr[,1],
  nTotComp=NA,
  n1=NA,
  n2=NA,
  delta=gr[,2],
  sd=gr[,3],
  es=NA,
  power=NA
)

datPowAnaem_sero$nTotComp<-floor((1-datPowAnaem_sero$pDropOut)*datPowAnaem_sero$nTot)
datPowAnaem_sero$n1<-round(datPowAnaem_sero$prev*datPowAnaem_sero$nTotComp)
datPowAnaem_sero$n2<-round((1-datPowAnaem_sero$prev)*datPowAnaem_sero$nTotComp)
datPowAnaem_sero$es<-round(digits=2,datPowAnaem_sero$delta/datPowAnaem_sero$sd)

for(j in 1:nrow(datPowAnaem_sero)){
  datPowAnaem_sero$power[j]<-pwr.t2n.test(n1=datPowAnaem_sero$n1[j],n2=datPowAnaem_sero$n2[j],d=datPowAnaem_sero$es[j],sig.level=datPowAnaem_sero$alpha[j])$power
}

datPowAnaem_sero$EsSd<-factor(paste(sep="","ES = ",datPowAnaem_sero$delta,", SD = ",datPowAnaem_sero$sd))

datPowAnaem_sero$EsSd
datPowAnaem_sero


pdf(file="~/work/PhD/VacciNTS/Sample size/SAINTS_protocol_sample_size/outputs/Vacc-iNTS_SampleSize_ananamia_seroconvert.pdf",width=16,height=9)
ggplot(data=datPowAnaem_sero,mapping=aes(x=nTot,y=power,group=EsSd, colour=EsSd)) +
  geom_point(size=3) +
  geom_line(lwd=1.5) +
  scale_color_manual(values=c("steelblue","greenyellow","salmon","orange"),name="effect size and std. dev.") +
  geom_hline(yintercept=0.8,lwd=1,lty=2,col="darkgrey") + 
  theme(legend.position = "right",text=element_text(size=16)) +
  xlab("Number of paired samples") +
  ylab("Power") +
  theme_bw() +
  ggtitle("Power curves; Seroconversion by anaemia status (prev 20%), drop-out (20%)")
dev.off()


## Sickle cell

# grid of sample size, effect size (proportion seroconverted in each group), standard deviation, prevalance of exposure
#  gr<-expand.grid(seq(1000,2500,by=100),c(0.1,0.2),c(0.4,0.6),c(1/24,0.074,0.08,0.099))
gr<-expand.grid(c(2000, 2225, 2500, 2750, 3000),c(0.2, 0.25),c(0.4,0.6),c(0.01))

datPowSickle_sero<-data.frame(
  alpha=0.05,
  prev=gr[,4],
  pDropOut=0.2,
  nTot=gr[,1],
  nTotComp=NA,
  n1=NA,
  n2=NA,
  delta=gr[,2],
  sd=gr[,3],
  es=NA,
  power=NA
)

datPowSickle_sero$nTotComp<-floor((1-datPowSickle_sero$pDropOut)*datPowSickle_sero$nTot)
datPowSickle_sero$n1<-round(datPowSickle_sero$prev*datPowSickle_sero$nTotComp)
datPowSickle_sero$n2<-round((1-datPowSickle_sero$prev)*datPowSickle_sero$nTotComp)
datPowSickle_sero$es<-round(digits=2,datPowSickle_sero$delta/datPowSickle_sero$sd)

for(j in 1:nrow(datPowSickle_sero)){
  datPowSickle_sero$power[j]<-pwr.t2n.test(n1=datPowSickle_sero$n1[j],n2=datPowSickle_sero$n2[j],d=datPowSickle_sero$es[j],sig.level=datPowSickle_sero$alpha[j])$power
}

datPowSickle_sero$EsSd<-factor(paste(sep="","ES = ",datPowSickle_sero$delta,", SD = ",datPowSickle_sero$sd))

datPowSickle_sero$EsSd
datPowSickle_sero


pdf(file="~/work/PhD/VacciNTS/Sample size/SAINTS_protocol_sample_size/outputs/Vacc-iNTS_SampleSize_sicklecell_seroconvert.pdf",width=16,height=9)
ggplot(data=datPowSickle_sero,mapping=aes(x=nTot,y=power,group=EsSd, colour=EsSd)) +
  geom_point(size=3) +
  geom_line(lwd=1.5) +
  scale_color_manual(values=c("steelblue","greenyellow","salmon","orange"),name="effect size and std. dev.") +
  geom_hline(yintercept=0.8,lwd=1,lty=2,col="darkgrey") + 
  theme(legend.position = "right",text=element_text(size=16)) +
  xlab("number of paired samples") +
  ylab("power") +
  theme_bw() +
  ggtitle("Power curves; Seroconversion by Sickle Cell status (prev 1%), drop-out (20%)")
dev.off()


#####################################################
########################
## Anaemia exposure
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


if(is.na(args[2]) | args[2]=="anaemia"){
  # simulate data
  simulateDataAnmSBA<-function(n,prevBin=0.2,sbaSD=0.4,sbaSlopeMonth=(-0.1),iccSD=0.03,effectSize=0.5,effectSizeSD=0.25,pDropOut=0.2,ageMin=8,ageMax=36){
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
    
    datSim$sba<-datSim$sba8m+datSim$u+(sbaSlopeMonth+datSim$v)*(datSim$ageSlopeEffectLimit-9)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
    datSim$sba[datSim$visit==2 & rbinom(n=n,size=1,prob=pDropOut)==1]<-NA
    
    return(datSim)
  }
  
  # detect a change in slope of SBA with age for i) malaria cases, ii) malaria vaccinated, iii) other binary covariates (linear mixed model approach)
  fitGEETestCoeff<-function(dat,alpha=0.05,coefName1="exposureexposed",coefName2=NULL){
    library(gee)
    #mod<-lmer(sba~agem9*exposure+(1|pid),REML=F,data=dat) # REML=F as we test fixed parameters; cannot have a second random effect even though this is how data was generated; would need more data points per participant...
    mod<-gee(sba~agem8*exposure,id=as.factor(pid),data=dat[order(dat$pid,dat$visit),],corstr="exchangeable",na.action=na.omit)
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
  
  
  # Anaemia; prevalence ~14%; try sample sizes 1000, 1500, 2000 and effect sizes 0.05, 0.1, 0.15
  gr<-expand.grid(c(2000, 2500, 3000),c(0.2, 0.225, 0.25, 0.275, 0.3),c(0.2))
  datPow<-data.frame(
    alpha=0.05,
    AnaemiaPrev=0.2,
    effectSize=gr[,2],
    n=gr[,1],
    pDropOut=gr[,3],
    powerShift=NA,
    powerSlope=NA
  )
  
  for(j in 1:nrow(datPow)){
    cat(paste(sep="","## effectSize = ",datPow$effectSize[j],", n = ",datPow$n[j],", pDropOut = ",datPow$pDropOut[j]," ##\n"))
    
    
    res<-foreach(i = 1:B) %dopar% {

      
      datSim<-simulateDataAnmSBA(n=datPow$n[j],prevBin=datPow$AnaemiaPrev[j],sbaSD=0.4,sbaSlopeMonth=(-0.1),iccSD=0.03,effectSize=datPow$effectSize[j],effectSizeSD=0.25,pDropOut=datPow$pDropOut[j],ageMin=8,ageMax=36)
      tmp<-fitGEETestCoeff(dat=datSim,coefName1="exposureexposed",coefName2="agem9:exposureexposed",alpha=datPow$alpha[j])
 
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

datPowAnem <- datPow
export(datPow, here("data", "datPowAnam.RData"))
datPowAnaem <- import(here("data", "datPowAnam.RData"))

# produce power curves
pdf(file="~/work/PhD/VacciNTS/Sample size/Vacc-iNTS_SampleSize_MalariaInfection_powerCurves_pDO0.2.pdf",width=16,height=9)
ggplot(data=datPowAnaem[datPowAnaem$pDropOut==0.2,],mapping=aes(x=effectSize,y=powerShift,group=as.factor(n),color=as.factor(n))) +
  geom_point(size=3) +
  geom_line(lwd=1.5) +
  geom_hline(yintercept=0.8,lwd=2,lty=2,col="darkgrey") + 
  scale_color_manual(values=colorRampPalette(brewer.pal(9,"Paired"))(length(levels(factor(datPow$n)))),name="number of paired samples") +
  ylab("power") +
  xlab("effect size (shift in average SBA)") +
  ggtitle("Power curves: Anaemia and SBA for different effect sizes and sample sizes, assuming 20% drop-out.") +
  theme(text=element_text(size=18)) +
  theme_bw()
dev.off()


datSim<-simulateDataAnmSBA(n=1000,prevBin=0.20,sbaSD=0.4,sbaSlopeMonth=(-0.1),iccSD=0.03,effectSize=0.1,effectSizeSD=0.25,pDropOut=0.2,ageMin=9,ageMax=36)
mod<-gee(sba~agem9*exposure,id=as.factor(pid),data=datSim[order(datSim$pid,datSim$visit),],corstr="exchangeable",na.action=na.omit)

g<-ggplot(data=datSim,mapping=aes(x=agem9,y=sba,col=exposure,group=pid)) +
  geom_point() + 
  geom_line() +
  scale_color_manual(values=c(c1,c2)) +
  geom_abline(lwd=2,intercept=coef(mod)["(Intercept)"]+coef(mod)["exposureexposed"],slope=coef(mod)["agem9"]+coef(mod)["agem9:exposureexposed"],col="orange") + 
  geom_abline(lwd=2,intercept=coef(mod)["(Intercept)"],slope=coef(mod)["agem9"],col="steelblue")
print(g)
dev.off()

##########################
## Malnutrition


if(is.na(args[2]) | args[2]=="malnutrition"){
  # simulate data
  simulateDataMalSBA<-function(n,prevBin=0.02,sbaSD=0.4,sbaSlopeMonth=(-0.1),iccSD=0.03,effectSize=0.5,effectSizeSD=0.25,pDropOut=0.2,ageMin=9,ageMax=36){
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
    
    datSim$sba<-datSim$sba9m+datSim$u+(sbaSlopeMonth+datSim$v)*(datSim$ageSlopeEffectLimit-9)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
    datSim$sba[datSim$visit==2 & rbinom(n=n,size=1,prob=pDropOut)==1]<-NA
    
    return(datSim)
  }
  
  # detect a change in slope of SBA with age for i) malaria cases, ii) malaria vaccinated, iii) other binary covariates (linear mixed model approach)
  fitGEETestCoeff<-function(dat,alpha=0.05,coefName1="exposureexposed",coefName2=NULL){
    library(gee)
    mod<-gee(sba~agem9*exposure,id=as.factor(pid),data=dat[order(dat$pid,dat$visit),],corstr="exchangeable",na.action=na.omit)
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
  
  
  # Malnutrition prevalence ~14%; try sample sizes 1000, 1500, 2000 and effect sizes 0.05, 0.1, 0.15
  gr<-expand.grid(c(2000, 2500, 3000),c(0.4, 0.5, 0.6, 0.7),c(0.2))
  datPow<-data.frame(
    alpha=0.05,
    MalnPrev=0.02,
    effectSize=gr[,2],
    n=gr[,1],
    pDropOut=gr[,3],
    powerShift=NA,
    powerSlope=NA
  )
  
  for(j in 1:nrow(datPow)){
    cat(paste(sep="","## effectSize = ",datPow$effectSize[j],", n = ",datPow$n[j],", pDropOut = ",datPow$pDropOut[j]," ##\n"))
    
   
    
    res<-foreach(i = 1:B) %dopar% {
  
      
      datSim<-simulateDataMalSBA(n=datPow$n[j],prevBin=datPow$MalnPrev[j],sbaSD=0.4,sbaSlopeMonth=(-0.1),iccSD=0.03,effectSize=datPow$effectSize[j],effectSizeSD=0.25,pDropOut=datPow$pDropOut[j],ageMin=9,ageMax=36)
      tmp<-fitGEETestCoeff(dat=datSim,coefName1="exposureexposed",coefName2="agem9:exposureexposed",alpha=datPow$alpha[j])
 
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

datPowMaln <- datPow
export(datPowMaln, here("data", "datPowMaln.RData"))


# produce power curves
pdf(file="~/work/PhD/VacciNTS/Sample size/Vacc-iNTS_SampleSize_MalariaInfection_powerCurves_pDO0.2.pdf",width=16,height=9)
ggplot(data=datPowMaln[datPowMaln$pDropOut==0.2,],mapping=aes(x=effectSize,y=powerShift,group=as.factor(n),color=as.factor(n))) +
  geom_point(size=3) +
  geom_line(lwd=1.5) +
  geom_hline(yintercept=0.8,lwd=2,lty=2,col="darkgrey") + 
  scale_color_manual(values=colorRampPalette(brewer.pal(9,"Paired"))(length(levels(factor(datPowMaln$n)))),name="number of paired samples") +
  ylab("power") +
  xlab("effect size (shift in average SBA)") +
  ggtitle("Power curves: Impact of Malnutrition on SBA, prev (2%), drop-out (20% )") +
  theme(text=element_text(size=18)) +
  theme_bw()
dev.off()



##########################
## Sickle Cell


if(is.na(args[2]) | args[2]=="sickle"){
  # simulate data
  simulateData<-function(n,prevBin=0.01,sbaSD=0.4,sbaSlopeMonth=(-0.1),iccSD=0.03,effectSize=0.5,effectSizeSD=0.25,pDropOut=0.2,ageMin=9,ageMax=36){
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
    
    datSim$sba<-datSim$sba9m+datSim$u+(sbaSlopeMonth+datSim$v)*(datSim$ageSlopeEffectLimit-9)+rnorm(n=nrow(datSim),mean=effectSize,sd=effectSizeSD)*ifelse(datSim$exposure=="exposed",1,0)
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
  
  
  # Sickle cell prevalence ~14%; try sample sizes 1000, 1500, 2000 and effect sizes 0.05, 0.1, 0.15
  gr<-expand.grid(c(2000, 2500, 3000),c(0.8,0.9,1),c(0.2))
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
    
   
    res<-foreach(i = 1:B) %dopar% {
     
      datSim<-simulateData(n=datPow$n[j],prevBin=datPow$sicklePrev[j],sbaSD=0.4,sbaSlopeMonth=(-0.1),iccSD=0.03,effectSize=datPow$effectSize[j],effectSizeSD=0.25,pDropOut=datPow$pDropOut[j],ageMin=9,ageMax=36)
      tmp<-fitGEETestCoeff(dat=datSim,coefName1="exposureexposed",coefName2="agem9:exposureexposed",alpha=datPow$alpha[j])
     
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

datPowSickle <- datPow
export(datPowSickle, here("data", "datPowSickle.RData"))



# produce power curves
pdf(file="~/work/PhD/VacciNTS/Sample size/Vacc-iNTS_SampleSize_MalariaInfection_powerCurves_pDO0.2.pdf",width=16,height=9)
ggplot(data=datPowSickle[datPowSickle$pDropOut==0.2,],mapping=aes(x=effectSize,y=powerShift,group=as.factor(n),color=as.factor(n))) +
  geom_point(size=3) +
  geom_line(lwd=1.5) +
  geom_hline(yintercept=0.8,lwd=2,lty=2,col="darkgrey") + 
  scale_color_manual(values=colorRampPalette(brewer.pal(9,"Paired"))(length(levels(factor(datPowMaln$n)))),name="number of paired samples") +
  ylab("power") +
  xlab("effect size (shift in average SBA)") +
  ggtitle("Power curves: Impact of Malnutrition on SBA, prev (2%), drop-out (20% )") +
  theme(text=element_text(size=18)) +
  theme_bw()
dev.off()



###############################################
## Difference in proportions of seroconversions
################################################

pwr.t2n.test(d=0.10,n1=1600,n2=400,alternative="greater")

pwr.2p2n.test(h=0.10,n1=1900,n2=100,sig.level=0.05,alternative="greater")



################################################################
## Sample size calculation - MPO with stool exposure ##
################################################################
library(pwr)

## Compare mean MPO between children exposed and unexposed to NTS ##
# d= effect size, n1 = number in sample 1, n2= number in sample 2, alternative = alternative hypothesis - two-sided, greater or less

pwr.t2n.test(d=0.25,n1=900,n2=100,alternative="greater")
pwr.2p2n.test(h=0.35,n1=1940,n2=60,sig.level=0.05,alternative="greater")
pwr.2p2n.test(h=0.30,n1=1900,n2=100,sig.level=0.05,alternative="greater")

pwr.f2.test