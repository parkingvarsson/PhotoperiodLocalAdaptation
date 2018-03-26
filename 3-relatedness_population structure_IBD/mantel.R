#! /usr/bin/Rscript --no-save --no-restore

library(vegan)
setwd("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/IBD/plot")
genetic=read.table("genetic-distances.txt",header=F)
geographic=read.table("geo-distances.txt",header=F)
IBD=read.table("distances.summary.txt",header=F)
names(IBD)=c("Distance","Fst")

mantel(genetic/(1-genetic),geographic,method="pearson",permutations=999)
###Mantel statistic r: 0.2098
###      Significance: 0.059

###Upper quantiles of permutations (null model):
###  90%   95% 97.5%   99%
###0.173 0.214 0.247 0.302
###Permutation: free
###Number of permutations: 999

mantel(xdis = genetic/(1 - genetic), ydis = geographic, method = "pearson",      permutations = 9999)

#Mantel statistic r: 0.2098
#      Significance: 0.047

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99%
#0.166 0.206 0.239 0.293
#Permutation: free
#Number of permutations: 9999


png("mantel.IBD.png",width = 6, height = 5, units = 'in', res=300)
par(mar=c(4,4.5,2,1))

plot(IBD$Distance,IBD$Fst/(1-IBD$Fst),xlab="Distance (km)",ylab=expression(F[ST]/1-F[ST]),pch=20,col="grey40",ylim=c(-0.009,0.003))
z=lm(IBD$Fst/(1-IBD$Fst)~IBD$Distance)
abline(z,col="red",lwd=2) # equivalent to abline(reg = z) or
#text(300,0.002,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
text(300,0.002,"r=0.2098")
p=0.047
ifelse(p<0.001,text(300,0.001,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(300,0.001,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
dev.off()





