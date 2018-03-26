#! /usr/bin/Rscript --no-save --no-restore

.libPaths("/pica/h1/jingwang/R/x86_64-redhat-linux-gnu-library/3.2")
####this script is used to make a plot of diversity and Tajima's D separately in four northern, two middle and six southern populations
library(RColorBrewer)
library(dplyr)
library(ggplot2)
colors <- brewer.pal(10,"Paired")

library(data.table)

setwd("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/summary_sig_test")

###To test the statistical significance for pi, Fst, H12, CLR between the chr10 region and genome-wide averages for northern, middle and southern populations separately

############################################################################################################################################
##1. H12
g_N_h12=data.frame()
g_M_h12=data.frame()
g_S_h12=data.frame()

###1.1. genome-wide
for (chr in 1:19) {
#N-h12
g_N_h12_chr=paste("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/maf005/pos/rm_het_hwe/H12/out/SwAsp.genome.h12.north.chr",chr,".200snps.txt",sep="")
g_N_h12_chr_data=fread(g_N_h12_chr,header=F)
names(g_N_h12_chr_data)=c("center","left","right","n_hap","H1","H2","H12","H2_H1")

if (chr == 3) {
g_N_h12_chr_data=g_N_h12_chr_data[!which(g_N_h12_chr_data$center > 22812000 & g_N_h12_chr_data$center < 22832000),]
}

if(chr == 10) {
g_N_h12_chr_data=g_N_h12_chr_data[!which(g_N_h12_chr_data$center > 16300000 & g_N_h12_chr_data$center < 16950000),]
}

g_N_h12=data.frame(rbind(g_N_h12,g_N_h12_chr_data))

#M-h12
g_M_h12_chr=paste("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/maf005/pos/rm_het_hwe/H12/out/SwAsp.genome.h12.mid.chr",chr,".200snps.txt",sep="")
g_M_h12_chr_data=fread(g_M_h12_chr,header=F)
names(g_M_h12_chr_data)=c("center","left","right","n_hap","H1","H2","H12","H2_H1")

if (chr == 3) {
g_M_h12_chr_data=g_M_h12_chr_data[!which(g_M_h12_chr_data$center > 22812000 & g_M_h12_chr_data$center < 22832000),]
}

if(chr == 10) {
g_M_h12_chr_data=g_M_h12_chr_data[!which(g_M_h12_chr_data$center > 16300000 & g_M_h12_chr_data$center < 16950000),]
}

g_M_h12=data.frame(rbind(g_M_h12,g_M_h12_chr_data))

##S-h12
g_S_h12_chr=paste("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/maf005/pos/rm_het_hwe/H12/out/SwAsp.genome.h12.south.chr",chr,".200snps.txt",sep="")
g_S_h12_chr_data=fread(g_S_h12_chr,header=F)
names(g_S_h12_chr_data)=c("center","left","right","n_hap","H1","H2","H12","H2_H1")

if (chr == 3) {
g_S_h12_chr_data=g_S_h12_chr_data[!which(g_S_h12_chr_data$center > 22812000 & g_S_h12_chr_data$center < 22832000),]
}

if(chr == 10) {
g_S_h12_chr_data=g_S_h12_chr_data[!which(g_S_h12_chr_data$center > 16300000 & g_S_h12_chr_data$center < 16950000),]
}

g_S_h12=data.frame(rbind(g_S_h12,g_S_h12_chr_data))

}


############################################################################################################################################

###1.2 selected chr10 region 
#chr10-N
chr10_N_h12=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/H12/out/SwAsp.chr10.7scaffolds.h12.north.200snps.txt",header=F)
names(chr10_N_h12)=c("center","left","right","n_hap","H1","H2","H12","H2_H1")
#chr10-M
chr10_M_h12=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/H12/out/SwAsp.chr10.7scaffolds.h12.mid.200snps.txt",header=F)
names(chr10_M_h12)=c("center","left","right","n_hap","H1","H2","H12","H2_H1")
#chr10-S
chr10_S_h12=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/H12/out/SwAsp.chr10.7scaffolds.h12.south.200snps.txt",header=F)
names(chr10_S_h12)=c("center","left","right","n_hap","H1","H2","H12","H2_H1")


############################################################################################################################################

###1.3 Summary

###1.3.1 genome-wide
g_N_h12$pop="North"
g_M_h12$pop="Mid"
g_S_h12$pop="South"

g_N_h12$group="genome-wide"
g_M_h12$group="genome-wide"
g_S_h12$group="genome-wide"

###1.3.2 selected region on chr10
chr10_N_h12$pop="North"
chr10_M_h12$pop="Mid"
chr10_S_h12$pop="South"

chr10_N_h12$group="chr10-700kb"
chr10_M_h12$group="chr10-700kb"
chr10_S_h12$group="chr10-700kb"

###1.3.3 summary from genome-wide and chr10

h12=rbind(g_N_h12,g_M_h12,g_S_h12,chr10_N_h12,chr10_M_h12,chr10_S_h12)
h12$pop=factor(h12$pop,levels=c("North","Mid","South"))
h12$group=factor(h12$group,levels=c("genome-wide","chr10-700kb"))

###1.3.4 Making the plot 

lvl1=c("chr10-700kb","genome-wide")
lvl2=c("North","Mid","South")
factor1=as.factor(c(rep("chr10-700kb",8948),rep("genome-wide",4227077),rep("chr10-700kb",8948),rep("genome-wide",4227077),rep("chr10-700kb",8948),rep("genome-wide",4227077)))
factor2=as.factor(c(rep("North",4236025),rep("Mid",4236025),rep("South",4236025)))
plotgrp=factor(paste(factor2,factor1),levels=c(sapply(lvl2,paste,lvl1)))

h12=c(chr10_N_h12$H12,g_N_h12$H12,chr10_M_h12$H12,g_M_h12$H12,chr10_S_h12$H12,g_S_h12$H12)

#######create the plot
pdf(file="H12.chr10_genome.genome.group.pdf",width=7,height=4)
par(las=1)
par(mar=c(3,10,1,1))
cols=c("firebrick1","dodgerblue2","dodgerblue4")
at_1=c(1:2,4:5,7:8)
at=c(1.5,4.5,7.5)
boxplot(h12~plotgrp,notch=TRUE,at=at_1,xaxt="n",cex.axis=1.5,frame.plot = FALSE,xlab="",ylab="",cex.lab=2.5,labs=2,col=c(cols[3],"grey20",cols[2],"grey50",cols[1],"grey80"),outline=FALSE,ylim=c(0,1))
axis(1,at=at,labels=lvl2,tick=T,cex.axis=2)
mtext("H12",,las=1,side=2,line=4,cex=2.5)
segments(1,0.9,1,0.92)
segments(1,0.92,2,0.92)
segments(2,0.9,2,0.92)
text(1.5,0.94,"***",cex=2)
segments(4,0.9,4,0.92)
segments(4,0.92,5,0.92)
segments(5,0.9,5,0.92)
text(4.5,0.94,"***",cex=2)
segments(7,0.9,7,0.92)
segments(7,0.92,8,0.92)
segments(8,0.9,8,0.92)
text(7.5,0.94,"n.s.",cex=2)
dev.off()


####H2/H1
h21=c(chr10_N_h12$H2_H1,g_N_h12$H2_H1,chr10_M_h12$H2_H1,g_M_h12$H2_H1,chr10_S_h12$H2_H1,g_S_h12$H2_H1)

wilcox.test(g_N_h12$H2_H1,chr10_N_h12$H2_H1)

wilcox.test(g_M_h12$H2_H1,chr10_M_h12$H2_H1)

wilcox.test(g_S_h12$H2_H1,chr10_S_h12$H2_H1)

#######create the plot
pdf(file="H2_H1.chr10_genome.genome.group.pdf",width=7,height=4)
par(las=1)
par(mar=c(3,10,1,1))
cols=c("firebrick1","dodgerblue2","dodgerblue4")
at_1=c(1:2,4:5,7:8)
at=c(1.5,4.5,7.5)
boxplot(h21~plotgrp,notch=TRUE,at=at_1,xaxt="n",cex.axis=1.5,frame.plot = FALSE,xlab="",ylab="",cex.lab=2.5,labs=2,col=c(cols[3],"grey20",cols[2],"grey50",cols[1],"grey80"),outline=FALSE,ylim=c(0,1.1))
axis(1,at=at,labels=lvl2,tick=T,cex.axis=2)
mtext("H2/H1",,las=1,side=2,line=3.5,cex=2.5)
segments(1,1.02,1,1,04)
segments(1,1.04,2,1.04)
segments(2,1.02,2,1.04)
text(1.5,1.06,"***",cex=2)
segments(4,1.02,4,1.04)
segments(4,1.04,5,1.04)
segments(5,1.02,5,1.04)
text(4.5,1.06,"***",cex=2)
segments(7,1.02,7,1.04)
segments(7,1.04,8,1.04)
segments(8,1.02,8,1.04)
text(7.5,1.06,"n.s.",cex=2)
dev.off()

###1.3.5 Summary of the statistics
#north
median(g_N_h12$H12,na.rm=T)
quantile(g_N_h12$H12,na.rm=T,c(0.025,0.975))
median(chr10_N_h12$H12,na.rm=T)
quantile(chr10_N_h12$H12,na.rm=T,c(0.025,0.975))

wilcox.test(g_N_h12$H12,chr10_N_h12$H12)

##mid
median(g_M_h12$H12,na.rm=T)
quantile(g_M_h12$H12,na.rm=T,c(0.025,0.975))
median(chr10_M_h12$H12,na.rm=T)
quantile(chr10_M_h12$H12,na.rm=T,c(0.025,0.975))

wilcox.test(g_M_h12$H12,chr10_M_h12$H12)

#south
median(g_S_h12$H12,na.rm=T)
quantile(g_S_h12$H12,na.rm=T,c(0.025,0.975))
median(chr10_S_h12$H12,na.rm=T)
quantile(chr10_S_h12$H12,na.rm=T,c(0.025,0.975))

wilcox.test(g_S_h12$H12,chr10_S_h12$H12)

################statistics for H2/H1

#north
median(g_N_h12$H2_H1,na.rm=T)
quantile(g_N_h12$H2_H1,na.rm=T,c(0.025,0.975))
median(chr10_N_h12$H2_H1,na.rm=T)
quantile(chr10_N_h12$H2_H1,na.rm=T,c(0.025,0.975))

wilcox.test(g_N_h12$H2_H1,chr10_N_h12$H2_H1)

##mid
median(g_M_h12$H2_H1,na.rm=T)
quantile(g_M_h12$H2_H1,na.rm=T,c(0.025,0.975))
median(chr10_M_h12$H2_H1,na.rm=T)
quantile(chr10_M_h12$H2_H1,na.rm=T,c(0.025,0.975))

wilcox.test(g_M_h12$H2_H1,chr10_M_h12$H2_H1)

#south
median(g_S_h12$H12,na.rm=T)
quantile(g_S_h12$H12,na.rm=T,c(0.025,0.975))
median(chr10_S_h12$H12,na.rm=T)
quantile(chr10_S_h12$H12,na.rm=T,c(0.025,0.975))

wilcox.test(g_S_h12$H12,chr10_S_h12$H12)




