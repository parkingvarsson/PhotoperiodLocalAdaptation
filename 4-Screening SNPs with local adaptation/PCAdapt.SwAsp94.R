#! /usr/bin/Rscript --no-save --no-restore

library(data.table)
library(pcadapt)
###qvalue
.libPaths("/pica/h1/jingwang/R/x86_64-redhat-linux-gnu-library/3.1")
library(qqman)
library(qvalue)

setwd("/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/PCAdapt/")

#data=read4pcadapt("/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/LEA/SwAsp.94samples.lfmm")
data=fread("/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/LEA/SwAsp.94samples.lfmm")

###1. Choice of K
k=pcadapt(data,K=10)

pdf("SwAsp94.pcadapt.pdf")
plot(k,option="screeplot")
dev.off()

###score plot
poplist=c(rep(1,6),rep(2,7),rep(3,9),rep(4,9),rep(5,8),rep(6,9),rep(7,10),rep(8,9),rep(9,3),5,rep(9,5),rep(10,8),rep(11,4),rep(12,6))
popsize=c(6,7,9,9,9,9,10,9,8,8,4,6)
pdf("SwAsp94.pcadapt.scores.pdf")
plot(k,option="scores",pop=poplist)
dev.off()


###calculating Weir & Cockerham's FST
fst=fstCalc(data,poplist,12,popsize)
write.table(fst,file="pcadapt.fst.txt",sep="\t",quote=F,row.names=F,col.names=T)

###write p-values to the table
###method-1: Using the correlation rho between the jth SNP and the kth principal component to detect the snps associated with local adaptation
x <- pcadapt(data,K=1)
qobj=qvalue(x$pvalues)
snp_pcadapt=read.table("/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/PCAdapt/snp",header=T)
snp_pcadapt$pvalue=x$pvalues
snp_pcadapt$qvalue=qobj$qvalues
loadings=x$loadings[,1]
snp_pcadapt$loadings=loadings
snp_pcadapt$chi2_stat=x$chi2_stat

write.table(snp_pcadapt,file="pcadapt.p_qvalues.K1.txt",sep="\t",quote=F,row.names=F,col.names=T)

###method-2: Using the communality method (correspond to the proportion of variance of a SNP that is explained by the first K PCs)

x_com=pcadapt(data,K=1,method="communality")

pdf("pcadapt.K1.communality.stat.pdf")
plot(x_com,option="stat.distribution")
dev.off()

qobj=qvalue(x_com$pvalues)
snp_pcadapt_com=read.table("/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/PCAdapt/snp",header=T)
snp_pcadapt_com$pvalue=x_com$pvalues
snp_pcadapt_com$qvalue=qobj$qvalues
loadings=x_com$loadings[,1]
snp_pcadapt_com$loadings=loadings
snp_pcadapt_com$chi2_stat=x_com$chi2_stat

write.table(snp_pcadapt_com,file="pcadapt.p_qvalues.K1.com.txt",sep="\t",quote=F,row.names=F,col.names=T)



#x2 <- pcadapt(data,K=2)
#write.table(x2$pvalues,file="pcadapt.pvalues.K2.txt",sep="\t",quote=F,row.names=F,col.names=T)



