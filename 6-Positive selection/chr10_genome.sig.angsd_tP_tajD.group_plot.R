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
##1. Diversity (10kb window, 5kb sliding)

###1.1. genome-wide
g_N_pi=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/ANGSD/SFS/northern/northern.thetas.summary.10kb_5kb.txt",header=T)
g_M_pi=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/ANGSD/SFS/middle/middle.thetas.summary.10kb_5kb.txt",header=T)
g_S_pi=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/ANGSD/SFS/southern/southern.thetas.summary.10kb_5kb.txt",header=T)

############################################################################################################################################

###1.2 selected chr10 region 
chr10_N_pi=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/diversity/ANGSD/chr10_region/pop/northern/northern.chr10.thetas.summary.10kb_5kb.txt",header=T)
chr10_M_pi=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/diversity/ANGSD/chr10_region/pop/middle/middle.chr10.thetas.summary.10kb_5kb.txt",header=T)
chr10_S_pi=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/diversity/ANGSD/chr10_region/pop/southern/southern.chr10.thetas.summary.10kb_5kb.txt",header=T)

############################################################################################################################################

###1.3 Summary
g_N_pi_new=filter(g_N_pi,!Chr %in% c("Potra000799","Potra000908","Potra000342","Potra001246","Potra004002","Potra003230","Potra000530"))
g_M_pi_new=filter(g_M_pi,!Chr %in% c("Potra000799","Potra000908","Potra000342","Potra001246","Potra004002","Potra003230","Potra000530"))
g_S_pi_new=filter(g_S_pi,!Chr %in% c("Potra000799","Potra000908","Potra000342","Potra001246","Potra004002","Potra003230","Potra000530"))

###1.3.1 genome-wide
g_N_pi_new$pop="North"
g_M_pi_new$pop="Mid"
g_S_pi_new$pop="South"

g_N_pi_new$group="genome-wide"
g_M_pi_new$group="genome-wide"
g_S_pi_new$group="genome-wide"

###1.3.2 selected region on chr10
chr10_N_pi$pop="North"
chr10_M_pi$pop="Mid"
chr10_S_pi$pop="South"

chr10_N_pi$group="chr10-700kb"
chr10_M_pi$group="chr10-700kb"
chr10_S_pi$group="chr10-700kb"

###1.3.3 summary from genome-wide and chr10

angsd_sum=rbind(g_N_pi_new,g_M_pi_new,g_S_pi_new,chr10_N_pi,chr10_M_pi,chr10_S_pi)
angsd_sum$pop=factor(angsd_sum$pop,levels=c("North","Mid","South"))
angsd_sum$group=factor(angsd_sum$group,levels=c("genome-wide","chr10-700kb"))

###1.3.4 Making the plot 
#diversity
lvl1=c("chr10-700kb","genome-wide")
lvl2=c("North","Mid","South")
factor1=as.factor(c(rep("chr10-700kb",117),rep("genome-wide",38003),rep("chr10-700kb",117),rep("genome-wide",38003),rep("chr10-700kb",117),rep("genome-wide",38001)))
factor2=as.factor(c(rep("North",38120),rep("Mid",38120),rep("South",38118)))
plotgrp=factor(paste(factor2,factor1),levels=c(sapply(lvl2,paste,lvl1)))

pi=c(chr10_N_pi$tP.norm,g_N_pi_new$tP.norm,chr10_M_pi$tP.norm,g_M_pi_new$tP.norm,chr10_S_pi$tP.norm,g_S_pi_new$tP.norm)

#######create the plot
pdf(file="pi.chr10_genome.genome.group.pdf",width=7,height=4)
par(las=1)
par(mar=c(3,10,1,1))
cols=c("firebrick1","dodgerblue2","dodgerblue4")
at_1=c(1:2,4:5,7:8)
at=c(1.5,4.5,7.5)
boxplot(pi~plotgrp,notch=TRUE,at=at_1,xaxt="n",cex.axis=1.5,frame.plot = FALSE,xlab="",ylab="",cex.lab=2.5,labs=2,col=c(cols[3],"grey20",cols[2],"grey50",cols[1],"grey80"),outline=FALSE,ylim=c(0,0.035))
axis(1,at=at,labels=lvl2,tick=T,cex.axis=2)
mtext(expression(pi),las=1,side=2,line=6,cex=2.5)
#mtext("***",1,line=0.2,at=c(1,4,7))
segments(1,0.032,1,0.033)
segments(1,0.033,1.95,0.033)
segments(1.95,0.032,1.95,0.033)
text(1.5,0.034,"***",cex=2)
segments(4,0.032,4,0.033)
segments(4,0.033,5,0.033)
segments(5,0.032,5,0.033)
text(4.5,0.034,"***",cex=2)
segments(7,0.032,7,0.033)
segments(7,0.033,8,0.033)
segments(8,0.032,8,0.033)
text(7.5,0.034,"***",cex=2)
dev.off()


###1.3.5 Determine the significance differences between different groups of measurements
#north
median(g_N_pi_new$tP.norm,na.rm=T)
quantile(g_N_pi_new$tP.norm,na.rm=T,c(0.025,0.975))
median(chr10_N_pi$tP.norm,na.rm=T)
quantile(chr10_N_pi$tP.norm,na.rm=T,c(0.025,0.975))

wilcox.test(g_N_pi_new$tP.norm,chr10_N_pi$tP.norm)
#mid
median(g_M_pi_new$tP.norm,na.rm=T)
quantile(g_M_pi_new$tP.norm,na.rm=T,c(0.025,0.975))
median(chr10_M_pi$tP.norm,na.rm=T)
quantile(chr10_M_pi$tP.norm,na.rm=T,c(0.025,0.975))

wilcox.test(g_M_pi_new$tP.norm,chr10_M_pi$tP.norm)
#south
median(g_S_pi_new$tP.norm,na.rm=T)
quantile(g_S_pi_new$tP.norm,na.rm=T,c(0.025,0.975))
median(chr10_S_pi$tP.norm,na.rm=T)
quantile(chr10_S_pi$tP.norm,na.rm=T,c(0.025,0.975))

wilcox.test(g_S_pi_new$tP.norm,chr10_S_pi$tP.norm)




