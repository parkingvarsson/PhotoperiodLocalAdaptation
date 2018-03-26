setwd("~/Dropbox/PaperIII/data/population_structure/IBD/")
#install.packages("vegan")
library(vegan)

IBD=read.table("distances.summary.txt")
names(IBD)=c("Distance","Fst")
mantel(IBD$Fst/(1-IBD$Fst),IBD$Distance,method="pearson",permutations=999)

#png(filename="SwAsp.94samples.IBD.png",width = 5, height = 4, units = 'in', res=300)
pdf(file="SwAsp.94samples.IBD.pdf",width=6,height=5)
par(mar=c(4,4.5,2,1))

plot(IBD$Distance,IBD$Fst/(1-IBD$Fst),xlab="Distance (km)",ylab=expression(F[ST]/1-F[ST]),pch=20,cex=1.5,col="grey20",ylim=c(-0.009,0.003))
z=lm(IBD$Fst/(1-IBD$Fst)~IBD$Distance)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$Gene_num,use = "complete"),2)
#r=round(summary(z)$adj.r.squared,3)
r=cor.test(IBD$Distance,IBD$Fst/(1-IBD$Fst))
#p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col="dark blue",lwd=2) # equivalent to abline(reg = z) or
#text(300,0.002,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
text(300,0.002,"r=0.2098")

###mantel test, permutation p-value
p=0.0472
ifelse(p<0.001,text(300,0.001,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(300,0.001,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
dev.off()
