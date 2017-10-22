#--some functions for plotting RiboTagger time series data
#--this file by RW on 17 October 2017
#--IAU workshop participants: please download this file, and use 'source(file="<this file name>") to load these functions into your workspace

calculate.relative.abudance.within.samples<-function(m)
{
 res=NULL
 mcs<-colSums(m)
 for(c in 1:ncol(m))
 {
  res<-cbind(res,m[,c]/mcs[c])
 }
 rownames(res)<-rownames(m)
 colnames(res)<-colnames(m)
 return(res)
}

plot.series<-function(xdata,expMat,annoMat,fileName="panel.pdf")
{
 #--some example plots from RiboTagger .xls output files (can also be used for 16S amp seq data with same format)
 #--this function does not assume that the row order of expMat and annoMat are the same; rows are selected by rowname
 pdf(file=fileName)
 for(c in (1:nrow(expMat)))
 {
  curtag<-rownames(expMat)[c]
  curmain<-paste(colnames(annoMat),as.vector(annoMat[curtag,]),sep="__")
  matplot(xdata,t(expMat),type="l",col="lightgrey",lty=1,las=1,main=curtag,xlab="Sampling day",ylab="Relative abundance",xaxt="n")
  rug(side=1,x=xdata,ticksize=0.015)
  axis(side=1,at=xdata,labels=as.character(xdata),cex.axis=0.75,las=3)
  lines(xdata,expMat[curtag,],col=1,lty=2)
  legend(min(xdata),max(expMat),curmain,ncol=1,bty="n",cex=0.8)
 }
 dev.off()
}

