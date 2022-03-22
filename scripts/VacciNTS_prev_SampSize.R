setwd("C:/Users/mhenrion.MLW/work/MLW_LSTM")
#source("scripts/standardSampleSizeCalcs.R")

x<-seq(0,1,by=0.01)
y<-rep(NA,length(x))
yl<-rep(NA,length(x))
yll<-rep(NA,length(x))
ylu<-rep(NA,length(x))
yu<-rep(NA,length(x))
yul<-rep(NA,length(x))
yuu<-rep(NA,length(x))
moe<-rep(NA,length(x))
moel<-rep(NA,length(x))
moeu<-rep(NA,length(x))
moeAlt1<-rep(NA,length(x))
moeAlt2<-rep(NA,length(x))

n<-200
nAlt1<-100
nAlt2<-400
nSim<-1e5

for(i in 1:length(x)){
    print(paste(sep="/",i,length(x)))
    tmp_y<-rep(NA,nSim)
    tmp_yu<-rep(NA,nSim)
    tmp_yl<-rep(NA,nSim)
    tmp_moeAlt1<-rep(NA,nSim)
    tmp_moeAlt2<-rep(NA,nSim)
    
    for(j in 1:nSim){
        samp<-rbinom(n=1,size=n,p=x[i])
        binEst<-binom.test(x=samp,n=n)
        tmp_y[j]<-binEst$estimate
        tmp_yl[j]<-binEst$conf.int[1]
        tmp_yu[j]<-binEst$conf.int[2]

        samp<-rbinom(n=1,size=nAlt1,p=x[i])
        binEst<-binom.test(x=samp,n=nAlt1)
        tmp_moeAlt1[j]<-(binEst$conf.int[2]-binEst$conf.int[1])/2
        samp<-rbinom(n=1,size=n,p=x[i])
        binEst<-binom.test(x=samp,n=nAlt2)
        tmp_moeAlt2[j]<-(binEst$conf.int[2]-binEst$conf.int[1])/2
    }
    
    tmp_moe<-(tmp_yu-tmp_yl)/2

    y[i]<-mean(tmp_y)
    yl[i]<-quantile(tmp_y,probs=0.025); yll[i]<-quantile(tmp_yl,probs=0.025); ylu[i]<-quantile(tmp_yl,probs=0.975)
    yu[i]<-quantile(tmp_y,probs=0.975); yul[i]<-quantile(tmp_yu,probs=0.025); yuu[i]<-quantile(tmp_yu,probs=0.975)
    moe[i]<-mean(tmp_moe); moel[i]<-quantile(tmp_moe,probs=0.025); moeu[i]<-quantile(tmp_moe,probs=0.975)
    moeAlt1[i]<-mean(tmp_moeAlt1)
    moeAlt2[i]<-mean(tmp_moeAlt2)
}

png(paste(sep="","GordonMelita_Vacc-iNTS/moeSampSize_n",n,"_nSim",nSim,".png"),width=16,height=9,res=450,units="in")
par(mfrow=c(1,2))

plot(x,y,xlab="prevalence",ylab="estimated prevalence",col="black",type="l",lwd=2,asp=1,main="n=200")
#tmpCol<-col2rgb("lightblue")
#polygon(x=c(x,x[length(x):1]),y=c(yll,ylu[length(x):1]),col=rgb(maxColorValue=255,red=tmpCol[1],green=tmpCol[2],blue=tmpCol[3],alpha=230),border=NA)
#tmpCol<-col2rgb("lightblue")
#polygon(x=c(x,x[length(x):1]),y=c(yul,yuu[length(x):1]),col=rgb(maxColorValue=255,red=tmpCol[1],green=tmpCol[2],blue=tmpCol[3],alpha=230),border=NA)
tmpCol<-col2rgb("steelblue")
polygon(x=c(x,x[length(x):1]),y=c(yl,yu[length(x):1]),col=rgb(maxColorValue=255,red=tmpCol[1],green=tmpCol[2],blue=tmpCol[3],alpha=200),border=NA)
lines(x,y,col="black",lwd=2)
#lines(lty=2,x,yu,col="steelblue",lwd=1.5)
#lines(lty=2,x,yl,col="steelblue",lwd=1.5)
#legend(x="bottomright",lty=c(1,2,1),lwd=c(2,1.5,10),col=c("black","steelblue",rgb(maxColorValue=255,red=tmpCol[1],green=tmpCol[2],blue=tmpCol[3],alpha=230)),legend=c("estimate","95% confidence interval","95% confidence interval of the lower and upper bounds of the CI"))
legend(x="bottomright",lty=c(1,1),lwd=c(2,10),col=c("black",rgb(maxColorValue=255,red=tmpCol[1],green=tmpCol[2],blue=tmpCol[3],alpha=200)),legend=c("estimate","95% confidence interval"))

plot(x,moe,col="black",type="l",lwd=2,main="n=200",xlab="prevalence",ylab="margin of error")
par(xpd=T)
tmpCol<-col2rgb("steelblue")
polygon(x=c(x,x[length(x):1]),y=c(moel,moeu[length(x):1]),col=rgb(maxColorValue=255,red=tmpCol[1],green=tmpCol[2],blue=tmpCol[3],alpha=200),border=NA)
lines(x,moe,col="black",lwd=2)
legend(x="topright",lty=c(1,1),col=c("black",rgb(maxColorValue=255,red=tmpCol[1],green=tmpCol[2],blue=tmpCol[3],alpha=200)),lwd=c(2,10),legend=c("margin of error","95% confidence interval"))
dev.off()


idx<-seq(1,length(x),by=5)
dat<-data.frame(x=x[idx],moe_n100=round(moeAlt1[idx],digits=2),moe_n200=round(moe[idx],digits=2),moe_n400=round(moeAlt2[idx],digits=2))
print(dat)
write.csv(dat,file=paste(sep="","GordonMelita_Vacc-iNTS/moeSampSize_n",n,"_nSim",nSim,".csv"))

save(list=ls(),file=paste(sep="","GordonMelita_Vacc-iNTS/moeSampSize_n",n,"_nSim",nSim,".RData"))
