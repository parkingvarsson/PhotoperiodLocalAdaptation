#! /usr/bin/Rscript --no-save --no-restore
library(geneplotter)
library(RColorBrewer)
pal <- colorRampPalette(c("light blue", "blue", "yellow"))

library(data.table)

setwd("/pica/v1/b2011141_nobackup/local_adaptation_paper/asp201_94/PCAdapt")
data=fread("pcadapt.p_qvalues.K1.fst.txt",header=T)



png(filename="pcadapt.fst.cor.png",width=6,height=6,units='in',res=300)
mat=matrix(c(1,0,2,3),2,byrow=TRUE)
layout=layout(mat,c(4,1),c(1,4))

par(mar=c(0.5,5,0,0))

squared_loadings=hist(data$chi2_stat,plot=FALSE,breaks=60)
barplot(squared_loadings$density,axes=FALSE,col="grey",space=0)

par(mar=c(5,5,0,0))
par(cex.lab=1.5)

colors_chi2_fst=densCols(data$chi2_stat,data$Fst,colramp=pal)
plot(data$chi2_stat,data$Fst,cex.lab=1.5,pch=19,cex.lab=1.1,col=colors_chi2_fst,cex=.5,xlab=expression(paste("Squared loadings ",rho[j1]^2)),ylab=expression(F[ST]))
outlier=data[which(data$qvalue<0.05),]
par(new=T)
points(outlier$chi2_stat,outlier$Fst,col="black",pch=19,cex=.5)
text(20,0.5,expression(paste(rho,"=0.230")^"***"))

par(mar=c(5,0.5,0,0))
fst=hist(data$Fst,plot=FALSE,breaks=60)
barplot(fst$density,axes=FALSE,col="grey",space=0,horiz=TRUE)
dev.off()


