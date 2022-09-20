##############################################################
## Sample size for seroconversion between risk factor groups ##
##############################################################

rm(list = ls())

library(pwr)
library(ggplot2)
library(rio)
library(here)
library(RColorBrewer)
#tinytex::install_tinytex()


################################
## 3. Salmonella stool exposure
################################

gr<-expand.grid(seq(1000,2500,by=100),c(0.1,0.15,0.2,0.25),c(0.2,0.3,0.4,0.5),c(0.4,0.06,0.08,0.099))

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
  geom_point(size=3) +
  geom_line(lwd=1.5) +
  scale_color_manual(values=colorRampPalette(brewer.pal(9,"Spectral"))(length(levels(datPowSeroConv$p1p2))),name="prop. seroconverted in each group") +
  geom_hline(yintercept=0.8,lwd=2,lty=2,col="darkgrey") + 
  geom_vline(xintercept=2000,lwd=2,lty=2,col="darkgrey") +
  theme(legend.position = "right",text=element_text(size=16)) +
  xlab("number of paired samples") +
  ylab("power") +
  theme_bw() +
  ggtitle("Power curves for different combinations of 3-month seroconversion rates. Prevalence of exposure = 8%, drop-out = 20%.")
dev.off()



