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
##1. SweepFinder2
g_N_sf2=data.frame()
g_M_sf2=data.frame()
g_S_sf2=data.frame()

###1.1. genome-wide
for (chr in 1:19) {
#N-sf2
if (chr <10) {
chrom=paste("0",chr,sep="")
#North
g_N_sf2_chr=paste("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/SweepFinder2/chromosome/SwAsp_94samples.filter.gt.beagle.ihs.north.tped.PseudoChr",chrom,".chr_tped.sf2.out",sep="")
g_N_sf2_chr_data=fread(g_N_sf2_chr,header=T)
g_N_sf2_chr_data$chr=chr
#Mid
g_M_sf2_chr=paste("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/SweepFinder2/chromosome/SwAsp_94samples.filter.gt.beagle.ihs.mid.tped.PseudoChr",chrom,".chr_tped.sf2.out",sep="")
g_M_sf2_chr_data=fread(g_M_sf2_chr,header=T)
g_M_sf2_chr_data$chr=chr
#South
g_S_sf2_chr=paste("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/SweepFinder2/chromosome/SwAsp_94samples.filter.gt.beagle.ihs.south.tped.PseudoChr",chrom,".chr_tped.sf2.out",sep="")
g_S_sf2_chr_data=fread(g_S_sf2_chr,header=T)
g_S_sf2_chr_data$chr=chr
} else {
chrom=chr
#North
g_N_sf2_chr=paste("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/SweepFinder2/chromosome/SwAsp_94samples.filter.gt.beagle.ihs.north.tped.PseudoChr",chrom,".chr_tped.sf2.out",sep="")
g_N_sf2_chr_data=fread(g_N_sf2_chr,header=T)
g_N_sf2_chr_data$chr=chr
#Mid
g_M_sf2_chr=paste("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/SweepFinder2/chromosome/SwAsp_94samples.filter.gt.beagle.ihs.mid.tped.PseudoChr",chrom,".chr_tped.sf2.out",sep="")
g_M_sf2_chr_data=fread(g_M_sf2_chr,header=T)
g_M_sf2_chr_data$chr=chr
#South
g_S_sf2_chr=paste("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/SweepFinder2/chromosome/SwAsp_94samples.filter.gt.beagle.ihs.south.tped.PseudoChr",chrom,".chr_tped.sf2.out",sep="")
g_S_sf2_chr_data=fread(g_S_sf2_chr,header=T)
g_S_sf2_chr_data$chr=chr
}

if (chr == 3) {
g_N_sf2_chr_data=g_N_sf2_chr_data[!which(g_N_sf2_chr_data$location > 22812000 & g_N_sf2_chr_data$location < 22832000),]
g_M_sf2_chr_data=g_M_sf2_chr_data[!which(g_M_sf2_chr_data$location > 22812000 & g_M_sf2_chr_data$location < 22832000),]
g_S_sf2_chr_data=g_S_sf2_chr_data[!which(g_S_sf2_chr_data$location > 22812000 & g_S_sf2_chr_data$location < 22832000),]
}

if(chr == 10) {
g_N_sf2_chr_data=g_N_sf2_chr_data[!which(g_N_sf2_chr_data$location > 16300000 & g_N_sf2_chr_data$location < 16950000),]
g_M_sf2_chr_data=g_M_sf2_chr_data[!which(g_M_sf2_chr_data$location > 16300000 & g_M_sf2_chr_data$location < 16950000),]
g_S_sf2_chr_data=g_S_sf2_chr_data[!which(g_S_sf2_chr_data$location > 16300000 & g_S_sf2_chr_data$location < 16950000),]
}

g_N_sf2=data.frame(rbind(g_N_sf2,g_N_sf2_chr_data))
g_M_sf2=data.frame(rbind(g_M_sf2,g_M_sf2_chr_data))
g_S_sf2=data.frame(rbind(g_S_sf2,g_S_sf2_chr_data))

}

############################################################################################################################################

###1.2 selected chr10 region 
#chr10-N
chr10_N_sf2=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/SweepFinder2/SwAsp_94samples.filter.gt.beagle.ihs.north.chr10.2kb.sf2.out",header=T)
chr10_N_sf2$chr=10
#chr10-M
chr10_M_sf2=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/SweepFinder2/SwAsp_94samples.filter.gt.beagle.ihs.mid.chr10.2kb.sf2.out",header=T)
chr10_M_sf2$chr=10
#chr10-S
chr10_S_sf2=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/SweepFinder2/SwAsp_94samples.filter.gt.beagle.ihs.south.chr10.2kb.sf2.out",header=T)
chr10_S_sf2$chr=10

############################################################################################################################################

###1.3 Summary

###1.3.1 genome-wide
g_N_sf2$pop="North"
g_M_sf2$pop="Mid"
g_S_sf2$pop="South"

g_N_sf2$group="genome-wide"
g_M_sf2$group="genome-wide"
g_S_sf2$group="genome-wide"

###1.3.2 selected region on chr10
chr10_N_sf2$pop="North"
chr10_M_sf2$pop="Mid"
chr10_S_sf2$pop="South"

chr10_N_sf2$group="chr10-700kb"
chr10_M_sf2$group="chr10-700kb"
chr10_S_sf2$group="chr10-700kb"

###1.3.3 summary from genome-wide and chr10

sf2=rbind(g_N_sf2,g_M_sf2,g_S_sf2,chr10_N_sf2,chr10_M_sf2,chr10_S_sf2)
sf2$pop=factor(sf2$pop,levels=c("North","Mid","South"))
sf2$group=factor(sf2$group,levels=c("genome-wide","chr10-700kb"))

###1.3.4 Making the plot 

lvl1=c("chr10-700kb","genome-wide")
lvl2=c("North","Mid","South")
factor1=as.factor(c(rep("chr10-700kb",337),rep("genome-wide",215668),rep("chr10-700kb",337),rep("genome-wide",215669),rep("chr10-700kb",337),rep("genome-wide",215669)))
factor2=as.factor(c(rep("North",216005),rep("Mid",216006),rep("South",216006)))
plotgrp=factor(paste(factor2,factor1),levels=c(sapply(lvl2,paste,lvl1)))

sf2=c(chr10_N_sf2$LR,g_N_sf2$LR,chr10_M_sf2$LR,g_M_sf2$LR,chr10_S_sf2$LR,g_S_sf2$LR)

pdf(file="CLR.chr10_genome.genome.group.pdf",width=7,height=4)
par(las=1)
par(mar=c(3,10,1,1))
cols=c("firebrick1","dodgerblue2","dodgerblue4")
at_1=c(1:2,4:5,7:8)
at=c(1.5,4.5,7.5)
boxplot(sf2~plotgrp,notch=TRUE,at=at_1,xaxt="n",cex.axis=1.5,frame.plot = FALSE,xlab="",ylab="",cex.lab=2.5,labs=2,col=c(cols[3],"grey20",cols[2],"grey50",cols[1],"grey80"),outline=FALSE,ylim=c(0,8))
axis(1,at=at,labels=lvl2,tick=T,cex.axis=2)
mtext("CLR",,las=1,side=2,line=4.5,cex=2.5)
segments(1,7,1,7.2)
segments(1,7.2,2,7.2)
segments(2,7,2,7.2)
text(1.5,7.4,"***",cex=2)
segments(4,7,4,7.2)
segments(4,7.2,5,7.2)
segments(5,7,5,7.2)
text(4.5,7.4,"n.s.",cex=2)
segments(7,7,7,7.2)
segments(7,7.2,8,7.2)
segments(8,7,8,7.2)
text(7.5,7.4,"n.s.",cex=2)
dev.off()



###1.3.5 summary
#north
median(g_N_sf2$LR,na.rm=T)
quantile(g_N_sf2$LR,na.rm=T,c(0.025,0.975))
median(chr10_N_sf2$LR,na.rm=T)
quantile(chr10_N_sf2$LR,na.rm=T,c(0.025,0.975))

wilcox.test(g_N_sf2$LR,chr10_N_sf2$LR)

#mid
median(g_M_sf2$LR,na.rm=T)
quantile(g_M_sf2$LR,na.rm=T,c(0.025,0.975))
median(chr10_M_sf2$LR,na.rm=T)
quantile(chr10_M_sf2$LR,na.rm=T,c(0.025,0.975))

wilcox.test(g_M_sf2$LR,chr10_M_sf2$LR)

#south
median(g_S_sf2$LR,na.rm=T)
quantile(g_S_sf2$LR,na.rm=T,c(0.025,0.975))
median(chr10_S_sf2$LR,na.rm=T)
quantile(chr10_S_sf2$LR,na.rm=T,c(0.025,0.975))

wilcox.test(g_S_sf2$LR,chr10_S_sf2$LR)



