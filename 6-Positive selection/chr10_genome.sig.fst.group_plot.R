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
##1. Fst

###1.1. genome-wide
g_NS_fst=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/maf005/pos/rm_het_hwe/vcftools/genome_NvsS.10_5kb.windowed.weir.fst",header=T)
g_NM_fst=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/maf005/pos/rm_het_hwe/vcftools/genome_NvsM.10_5kb.windowed.weir.fst",header=T)
g_SM_fst=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/maf005/pos/rm_het_hwe/vcftools/genome_SvsM.10_5kb.windowed.weir.fst",header=T)
pseudo_coord=read.table("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/maf005/pos/rm_het_hwe/vcftools/pseudo_coord.7scaffolds.txt",header=F)

############################################################################################################################################

###1.2 selected chr10 region 
chr10_NS_fst=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/vcftools/chr10_NvsS.10_5kb.windowed.weir.fst",header=T)
chr10_NM_fst=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/vcftools/chr10_NvsM.10_5kb.windowed.weir.fst",header=T)
chr10_SM_fst=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/vcftools/chr10_SvsM.10_5kb.windowed.weir.fst",header=T)

############################################################################################################################################

###1.3 Summary

g_NS_fst_new=g_NS_fst[!which((g_NS_fst$CHROM == "PseudoChr03" & g_NS_fst$BIN_START >=22812000 & g_NS_fst$BIN_END <=22832000) | (g_NS_fst$CHROM == "PseudoChr10" & g_NS_fst$BIN_START >= 16300000 & g_NS_fst$BIN_END <= 16950000)),]
g_NM_fst_new=g_NM_fst[!which((g_NM_fst$CHROM == "PseudoChr03" & g_NM_fst$BIN_START >=22812000 & g_NM_fst$BIN_END <=22832000) | (g_NM_fst$CHROM == "PseudoChr10" & g_NM_fst$BIN_START >= 16300000 & g_NM_fst$BIN_END <= 16950000)),]
g_SM_fst_new=g_SM_fst[!which((g_SM_fst$CHROM == "PseudoChr03" & g_SM_fst$BIN_START >=22812000 & g_SM_fst$BIN_END <=22832000) | (g_SM_fst$CHROM == "PseudoChr10" & g_SM_fst$BIN_START >= 16300000 & g_SM_fst$BIN_END <= 16950000)),]


###1.3.1 genome-wide
g_NS_fst_new$pop="N vs. S"
g_NM_fst_new$pop="N vs. M"
g_SM_fst_new$pop="S vs. M"

g_NS_fst_new$group="genome-wide"
g_NM_fst_new$group="genome-wide"
g_SM_fst_new$group="genome-wide"

###1.3.2 selected region on chr10
chr10_NS_fst$pop="N vs. S"
chr10_NM_fst$pop="N vs. M"
chr10_SM_fst$pop="S vs. M"

chr10_NS_fst$group="chr10-700kb"
chr10_NM_fst$group="chr10-700kb"
chr10_SM_fst$group="chr10-700kb"

###1.3.3 summary from genome-wide and chr10

fst=rbind(g_NS_fst_new,g_NM_fst_new,g_SM_fst_new,chr10_NS_fst,chr10_NM_fst,chr10_SM_fst)
fst$pop=factor(fst$pop,levels=c("N vs. S","N vs. M","S vs. M"))
fst$group=factor(fst$group,levels=c("genome-wide","chr10-700kb"))

fst_new=filter(fst,N_VARIANTS>10)

##############
g_NS_fst_new=filter(g_NS_fst_new,N_VARIANTS>10)   ##69190
g_NM_fst_new=filter(g_NM_fst_new,N_VARIANTS>10)   ##69185
g_SM_fst_new=filter(g_SM_fst_new,N_VARIANTS>10)   ##69190
chr10_NS_fst=filter(chr10_NS_fst,N_VARIANTS>10)   ##134
chr10_NM_fst=filter(chr10_NM_fst,N_VARIANTS>10)   ##134
chr10_SM_fst=filter(chr10_SM_fst,N_VARIANTS>10)   ##134

###1.3.4 Making the plot 
#Fst
####create groups
lvl1=c("chr10-700kb","genome-wide")
lvl2=c("N vs. S","N vs. M","S vs. M")
factor1=as.factor(c(rep("chr10-700kb",134),rep("genome-wide",69190),rep("chr10-700kb",134),rep("genome-wide",69185),rep("chr10-700kb",134),rep("genome-wide",69190)))
factor2=as.factor(c(rep("N vs. S",69324),rep("N vs. M",69319),rep("S vs. M",69324)))
plotgrp=factor(paste(factor2,factor1),levels=c(sapply(lvl2,paste,lvl1)))

fst=c(chr10_NS_fst$WEIGHTED_FST,g_NS_fst_new$WEIGHTED_FST,chr10_NM_fst$WEIGHTED_FST,g_NM_fst_new$WEIGHTED_FST,chr10_SM_fst$WEIGHTED_FST,g_SM_fst_new$WEIGHTED_FST)

####create the plot
pdf(file="fst.chr10_genome.genome.group.pdf",width=7,height=4)
par(las=1)
par(mar=c(3,10,1,1))
cols=c("firebrick1","dodgerblue2","dodgerblue4")
at_1=c(1:2,4:5,7:8)
at=c(1.5,4.5,7.5)
boxplot(fst~plotgrp,notch=TRUE,at=at_1,xaxt="n",cex.axis=1.5,frame.plot = FALSE,xlab="",ylab="",cex.lab=2.5,labs=2,col=c(cols[3],"grey20",cols[2],"grey50",cols[1],"grey80"),outline=FALSE,ylim=c(-0.05,0.5))
axis(1,at=at,labels=lvl2,tick=T,cex.axis=2)
mtext(expression(F[ST]),las=1,side=2,line=4,cex=2.5)
segments(1,0.45,1,0.46)
segments(1,0.46,2,0.46)
segments(2,0.45,2,0.46)
text(1.5,0.47,"***",cex=2)
segments(4,0.45,4,0.46)
segments(4,0.46,5,0.46)
segments(5,0.45,5,0.46)
text(4.5,0.47,"***",cex=2)
segments(7,0.45,7,0.46)
segments(7,0.46,8,0.46)
segments(8,0.45,8,0.46)
text(7.5,0.47,"***",cex=2)
dev.off()

###1.3.5 Summary
#N vs. S
median(g_NS_fst_new$WEIGHTED_FST,na.rm=T)
quantile(g_NS_fst_new$WEIGHTED_FST,na.rm=T,c(0.025,0.975))
median(chr10_NS_fst$WEIGHTED_FST,na.rm=T)
quantile(chr10_NS_fst$WEIGHTED_FST,na.rm=T,c(0.025,0.975))

wilcox.test(g_NS_fst_new$WEIGHTED_FST,chr10_NS_fst$WEIGHTED_FST)

#N vs. M
median(g_NM_fst_new$WEIGHTED_FST,na.rm=T)
quantile(g_NM_fst_new$WEIGHTED_FST,na.rm=T,c(0.025,0.975))
median(chr10_NM_fst$WEIGHTED_FST,na.rm=T)
quantile(chr10_NM_fst$WEIGHTED_FST,na.rm=T,c(0.025,0.975))

wilcox.test(g_NM_fst_new$WEIGHTED_FST,chr10_NM_fst$WEIGHTED_FST)

#S vs. M
median(g_SM_fst_new$WEIGHTED_FST,na.rm=T)
quantile(g_SM_fst_new$WEIGHTED_FST,na.rm=T,c(0.025,0.975))
median(chr10_SM_fst$WEIGHTED_FST,na.rm=T)
quantile(chr10_SM_fst$WEIGHTED_FST,na.rm=T,c(0.025,0.975))

wilcox.test(g_SM_fst_new$WEIGHTED_FST,chr10_SM_fst$WEIGHTED_FST)






