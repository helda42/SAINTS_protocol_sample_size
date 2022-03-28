##############################################################
## Sample size for seroconversion between risk factor groups ##
##############################################################

rm(list = ls())

library(pwr)
library(ggplot2)
#tinytex::install_tinytex()


## Malnutrition

# grid of sample size, effect size (proportion seroconverted in each group), standard deviation, prevalance of exposure
#  gr<-expand.grid(seq(1000,2500,by=100),c(0.1,0.2),c(0.4,0.6),c(1/24,0.074,0.08,0.099))
gr<-expand.grid(c(1500, 1750, 2000),c(0.1, 0.2),c(0.4,0.6),c(0.03))

datPowSalm<-data.frame(
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

datPowSalm$nTotComp<-floor((1-datPowSalm$pDropOut)*datPowSalm$nTot)
datPowSalm$n1<-round(datPowSalm$prev*datPowSalm$nTotComp)
datPowSalm$n2<-round((1-datPowSalm$prev)*datPowSalm$nTotComp)
datPowSalm$es<-round(digits=2,datPowSalm$delta/datPowSalm$sd)

for(j in 1:nrow(datPowSalm)){
  datPowSalm$power[j]<-pwr.t2n.test(n1=datPowSalm$n1[j],n2=datPowSalm$n2[j],d=datPowSalm$es[j],sig.level=datPowSalm$alpha[j])$power
}

datPowSalm$EsSd<-factor(paste(sep="","ES = ",datPowSalm$delta,", SD = ",datPowSalm$sd))

datPowSalm$EsSd
datPowSalm


pdf(file="~/work/PhD/VacciNTS/Sample size/SAINTS_protocol_sample_size/outputs/Vacc-iNTS_SampleSize_malnutrition_seroconvert.pdf",width=16,height=9)
ggplot(data=datPowSalm,mapping=aes(x=nTot,y=power,group=EsSd, colour=EsSd)) +
  geom_point(size=4) +
  geom_line(lwd=2) +
  scale_color_manual(values=c("steelblue","greenyellow","salmon","orange"),name="effect size and std. dev.") +
  geom_hline(yintercept=0.8,lwd=2,lty=2,col="darkgrey") + 
  theme(legend.position = "right",text=element_text(size=16)) +
  xlab("number of paired samples") +
  ylab("power") +
  ggtitle("Power curves by effect size. Prevalence of Malnutrition = 0.03, drop-out = 20%.") +
  theme_bw()
dev.off()


## Malaria

# grid of sample size, effect size (proportion seroconverted in each group), standard deviation, prevalance of exposure
#  gr<-expand.grid(seq(1000,2500,by=100),c(0.1,0.2),c(0.4,0.6),c(1/24,0.074,0.08,0.099))
gr<-expand.grid(c(1750, 2000, 2225, 2500),c(0.05, 0.1),c(0.4,0.6),c(0.16))

datPowSalm<-data.frame(
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

datPowSalm$nTotComp<-floor((1-datPowSalm$pDropOut)*datPowSalm$nTot)
datPowSalm$n1<-round(datPowSalm$prev*datPowSalm$nTotComp)
datPowSalm$n2<-round((1-datPowSalm$prev)*datPowSalm$nTotComp)
datPowSalm$es<-round(digits=2,datPowSalm$delta/datPowSalm$sd)

for(j in 1:nrow(datPowSalm)){
  datPowSalm$power[j]<-pwr.t2n.test(n1=datPowSalm$n1[j],n2=datPowSalm$n2[j],d=datPowSalm$es[j],sig.level=datPowSalm$alpha[j])$power
}

datPowSalm$EsSd<-factor(paste(sep="","ES = ",datPowSalm$delta,", SD = ",datPowSalm$sd))

datPowSalm$EsSd
datPowSalm


pdf(file="~/work/PhD/VacciNTS/Sample size/SAINTS_protocol_sample_size/outputs/Vacc-iNTS_SampleSize_malaria_seroconvert.pdf",width=16,height=9)
ggplot(data=datPowSalm,mapping=aes(x=nTot,y=power,group=EsSd, colour=EsSd)) +
  geom_point(size=4) +
  geom_line(lwd=2) +
  scale_color_manual(values=c("steelblue","greenyellow","salmon","orange"),name="effect size and std. dev.") +
  geom_hline(yintercept=0.8,lwd=2,lty=2,col="darkgrey") + 
  theme(legend.position = "right",text=element_text(size=16)) +
  xlab("number of paired samples") +
  ylab("power") +
  theme_bw() +
  ggtitle("Power curves by effect size. Prevalence of Malaria = 0.16, drop-out = 20%.")
dev.off()



## Anaemia

# grid of sample size, effect size (proportion seroconverted in each group), standard deviation, prevalance of exposure
#  gr<-expand.grid(seq(1000,2500,by=100),c(0.1,0.2),c(0.4,0.6),c(1/24,0.074,0.08,0.099))
gr<-expand.grid(c(1750, 2000, 2225, 2500),c(0.05, 0.1),c(0.4,0.6),c(0.34))

datPowSalm<-data.frame(
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

datPowSalm$nTotComp<-floor((1-datPowSalm$pDropOut)*datPowSalm$nTot)
datPowSalm$n1<-round(datPowSalm$prev*datPowSalm$nTotComp)
datPowSalm$n2<-round((1-datPowSalm$prev)*datPowSalm$nTotComp)
datPowSalm$es<-round(digits=2,datPowSalm$delta/datPowSalm$sd)

for(j in 1:nrow(datPowSalm)){
  datPowSalm$power[j]<-pwr.t2n.test(n1=datPowSalm$n1[j],n2=datPowSalm$n2[j],d=datPowSalm$es[j],sig.level=datPowSalm$alpha[j])$power
}

datPowSalm$EsSd<-factor(paste(sep="","ES = ",datPowSalm$delta,", SD = ",datPowSalm$sd))

datPowSalm$EsSd
datPowSalm


pdf(file="~/work/PhD/VacciNTS/Sample size/SAINTS_protocol_sample_size/outputs/Vacc-iNTS_SampleSize_ananamia_seroconvert.pdf",width=16,height=9)
ggplot(data=datPowSalm,mapping=aes(x=nTot,y=power,group=EsSd, colour=EsSd)) +
  geom_point(size=4) +
  geom_line(lwd=2) +
  scale_color_manual(values=c("steelblue","greenyellow","salmon","orange"),name="effect size and std. dev.") +
  geom_hline(yintercept=0.8,lwd=2,lty=2,col="darkgrey") + 
  theme(legend.position = "right",text=element_text(size=16)) +
  xlab("number of paired samples") +
  ylab("power") +
  theme_bw() +
  ggtitle("Power curves by effect size. Prevalence of Anaemia = 0.34, drop-out = 20%.")
dev.off()


## Sickle cell

# grid of sample size, effect size (proportion seroconverted in each group), standard deviation, prevalance of exposure
#  gr<-expand.grid(seq(1000,2500,by=100),c(0.1,0.2),c(0.4,0.6),c(1/24,0.074,0.08,0.099))
gr<-expand.grid(c(1750, 2000, 2225, 2500),c(0.1, 0.2),c(0.4,0.6),c(0.02))

datPowSalm<-data.frame(
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

datPowSalm$nTotComp<-floor((1-datPowSalm$pDropOut)*datPowSalm$nTot)
datPowSalm$n1<-round(datPowSalm$prev*datPowSalm$nTotComp)
datPowSalm$n2<-round((1-datPowSalm$prev)*datPowSalm$nTotComp)
datPowSalm$es<-round(digits=2,datPowSalm$delta/datPowSalm$sd)

for(j in 1:nrow(datPowSalm)){
  datPowSalm$power[j]<-pwr.t2n.test(n1=datPowSalm$n1[j],n2=datPowSalm$n2[j],d=datPowSalm$es[j],sig.level=datPowSalm$alpha[j])$power
}

datPowSalm$EsSd<-factor(paste(sep="","ES = ",datPowSalm$delta,", SD = ",datPowSalm$sd))

datPowSalm$EsSd
datPowSalm


pdf(file="~/work/PhD/VacciNTS/Sample size/SAINTS_protocol_sample_size/outputs/Vacc-iNTS_SampleSize_sicklecell_seroconvert.pdf",width=16,height=9)
ggplot(data=datPowSalm,mapping=aes(x=nTot,y=power,group=EsSd, colour=EsSd)) +
  geom_point(size=4) +
  geom_line(lwd=2) +
  scale_color_manual(values=c("steelblue","greenyellow","salmon","orange"),name="effect size and std. dev.") +
  geom_hline(yintercept=0.8,lwd=2,lty=2,col="darkgrey") + 
  theme(legend.position = "right",text=element_text(size=16)) +
  xlab("number of paired samples") +
  ylab("power") +
  theme_bw() +
  ggtitle("Power curves by effect size. Prevalence of Sickle Cell = 0.02, drop-out = 20%.")
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