setwd("~/Dropbox/PaperIII/data/local_adaptation_test/LFMM/env_data")
library("lattice")
#install.packages("corrplot")
library(corrplot)

env=read.table("SwAsp.Clone.Origin.data.txt",header=T)

env_cor=cor(env[,c(c(4:16),c(18:43))],method="spearman")
levelplot(env_cor)
corrplot(env_cor,method="circle")
par(mfrow=c(1,1))

pdf("env_corr.pdf")
corrplot(env_cor,type="upper",method="ellipse",tl.col="black", cl.cex = 1,tl.cex=0.6,tl.srt=45,order="hclust")
dev.off()

#use prcomp() and use center=TRUE to set the means of all variables equal to zero, and the scale=T, to set the standard deviations of all variables to be equal to 1.
pca=prcomp(env[,c(c(4:16),c(18:43))],scale=TRUE,center=TRUE)
sd=pca$sdev
summary_pca=summary(pca)

#write.table(as.data.frame(loadings[,1:4]),file=paste("environmental_pca.loadings.txt",sep="\t",quote=F,row.names=T,col.names=T))

#loadings
loadings=pca$rotation

#scores of each PC
scores=pca$x
scores_PC1_3=as.data.frame(scores[,c(1:3)])
scores_PC1_3$ind=env[,1]
scores_PC1_3$pop[scores_PC1_3$ind<11]="Pop1"
scores_PC1_3$pop[scores_PC1_3$ind>10&scores_PC1_3$ind<21]="Pop2"
scores_PC1_3$pop[scores_PC1_3$ind>20&scores_PC1_3$ind<31]="Pop3"
scores_PC1_3$pop[scores_PC1_3$ind>30&scores_PC1_3$ind<41]="Pop4"
scores_PC1_3$pop[scores_PC1_3$ind>40&scores_PC1_3$ind<51]="Pop5"
scores_PC1_3$pop[scores_PC1_3$ind>50&scores_PC1_3$ind<61]="Pop6"
scores_PC1_3$pop[scores_PC1_3$ind>60&scores_PC1_3$ind<71]="Pop7"
scores_PC1_3$pop[scores_PC1_3$ind>70&scores_PC1_3$ind<81]="Pop8"
scores_PC1_3$pop[scores_PC1_3$ind>80&scores_PC1_3$ind<91]="Pop9"
scores_PC1_3$pop[scores_PC1_3$ind>90&scores_PC1_3$ind<101]="Pop10"
scores_PC1_3$pop[scores_PC1_3$ind>100&scores_PC1_3$ind<111]="Pop11"
scores_PC1_3$pop[scores_PC1_3$ind>110&scores_PC1_3$ind<117]="Pop12"

#make the scores plot
col.rainbow <- rainbow(12)
palette(col.rainbow)
pop=as.factor(scores_PC1_3$pop)

#write the table for output
write.table(scores_PC1_3,file="environmental_pca.scores.txt",sep="\t",quote=F,row.names=F,col.names=T)
#write.table(scores_PC1_3,file="environmental_pca.scores_2.txt",sep="|",quote=F,row.names=F,col.names=T)

pdf("PC1.pdf")
plot(scores_PC1_3$ind,scores_PC1_3$PC1,col=pop,xlab="SwAsp samples",ylab="PC1")
dev.off()
pdf("PC1.latitude.pdf")
plot(env$Tree.latitude,scores_PC1_3$PC1,col=pop,xlab="Latitude",ylab="PC1")
dev.off()

pdf("PC2.pdf")
plot(scores_PC1_3$ind,scores_PC1_3$PC2,col=pop,xlab="SwAsp samples",ylab="PC2")
dev.off()
pdf("PC2,longitude.pdf")
plot(env$Tree.longitude,scores_PC1_3$PC2,col=pop,xlab="Longitude",ylab="PC2")
dev.off()

pdf("PC3.pdf")
plot(scores_PC1_3$ind,scores_PC1_3$PC3,col=pop,xlab="SwAsp samples",ylab="PC3")
dev.off()
pdf("PC3.elevation.pdf")
plot(env$Elevation,scores_PC1_3$PC3,col=pop,xlab="Elevation",ylab="PC3")
dev.off()


 #Variance and cutoff 
var=sd^2
var.percent=var/sum(var)*100
png("pca.percent.png",width = 5, height =4, units = 'in', res=300)
par(mar=c(4,5,1,1))
barplot(var.percent,xlab="Environmental PCs",ylab="Percent Variance",names.arg=1:length(var.percent),las=1,ylim=c(0,max(var.percent)),col="grey")
#abline(h=1/ncol(env)*100,col="red")
dev.off()

#cutoff for important loadings
sqrt(1/ncol(env))
loadings

#ilustrate the distribution of the samples in each ordination space
dev.new(height=7,width=7)

pdf("biplot_pc1_2.pdf")
biplot(scores[,1:2],loadings[,1:2],cex=0.7)
dev.off()

pdf("biplot_pc3_4.pdf")
biplot(scores[,3:4],loadings[,3:4],cex=0.5)
dev.off()




#using princomp()
pca1=princomp(env[,c(c(4:16),c(18:43))])
summary(pca1)
#sqrt of eigenvalues

biplot(pca1)
screeplot(pca1)
pca1$loadings
pca1$scores
loadings(pca1)
varimax(pca1$loadings[,1:3])

#use princomp with corr=TRUE
pca2=princomp(env[,c(c(4:16),c(18:43))],cor=TRUE)
summary(pca2)
biplot(pca2)
pca$sdev 
#loadings
unclass(pca2$loadings)
#PCs
pca2$scores

#4. Use PCA() in package "FactoMineR"
#install.packages("FactoMineR")
library(FactoMineR)
pca3=PCA(env[,c(c(4:16),c(18:43))],graph=FALSE)
#matrix with eigenvalues
pca3$eig
#correlations between variables and PCs
pca3$var$coord
#PCs
pca3$ind$coord


