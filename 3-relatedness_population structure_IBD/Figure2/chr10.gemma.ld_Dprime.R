#! /usr/bin/Rscript --no-save --no-restore
.libPaths("/pica/h1/jingwang/R/x86_64-redhat-linux-gnu-library/3.2")

library(data.table)
library("RColorBrewer")
source("/proj/b2011141/pipeline/R/local_adaptation_paper/manhattan/20171015/legend.col.R")

####################################################################
###Main aim: to create the GEMMA plot specifically for the Chr10 selected region, and also include the LD pattern (with the most significant SNP: Potra001245:25256) in it
setwd("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/sliding_three_methods")

chr10=fread("chr10.7scaffolds.region.three.txt",header=T)
ld<-fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/ld_r2_pos/SwAsp.chr10.7scaffolds.recode.Potri.rm_het.gt.beagle.hap.Dprimer.ld",header=T)

##Select LD for most associated SNP (SNP:Potra001246:25256)
ld[ld$POS1== 396298 |ld$POS2==396298,]->ldDp
ldDp$POS<-ifelse(ldDp$POS1==396298,ldDp$POS2,ldDp$POS1)
names(ldDp)[5]<-"R2"
ldDp<-rbind(ldDp,list("Chr10",396298,396298,188,1,1,1,396298))
ldDp<-setorder(ldDp,POS)
chr10$R2<-ldDp$R2
chr10$Dprime<-ldDp$Dprime

##Mark the two SNPs which were inferred to be causal SNP by CAVIAR, Potra001246:25256 and Potra001246:43095
caviar_snp=chr10[which(chr10$Physical=="378459" |chr10$Physical=="396298"),]

##Make the GWAS plot
pdf("GWAS_chr10_Dprime.pdf",width=12,height=6)
par(mgp=c(2.25,1,0))
par(mar=c(5,5,2,2))
colfunc <- colorRampPalette(c("dodgerblue4","darkorange","firebrick1"))
cols<-colfunc(10)[ceiling(abs(chr10$Dprime)*10)]
plot(-log10(GEMMA_p_lrt)~ Physical,data=chr10,yaxt="n",xaxt="n",bty="n",ylab=expression(-log[10](p)),xlab="Position",type="n",cex.axis=2,cex.lab=2)
points(-log10(GEMMA_p_lrt)~Physical,data=chr10,col=cols,pch=20)
points(-log10(GEMMA_p_lrt)~Physical,data=caviar_snp,col="black",pch=1,cex=1.5,lwd=2)
axis(1,lwd=2,at=seq(1e5*floor(min(chr10$Physical)/1e5),1e5*ceiling(max(chr10$Physical)/1e5),0.1*1e6),labels=seq(0,700,100),cex.axis=1.5)
axis(2,las=1,lwd=2,cex.axis=1.5)
legend.col(col = colfunc(100), lev = chr10$Dprime)
dev.off()

pdf("GWAS_chr10_R2.pdf",width=12,height=6)
par(mgp=c(2.25,1,0))
colfunc <- colorRampPalette(c("dodgerblue4","darkorange","firebrick1"))
cols<-colfunc(10)[ceiling(abs(chr10$R2)*10)]
plot(-log10(GEMMA_p_lrt)~ Physical,data=chr10,yaxt="n",xaxt="n",bty="n",ylab=expression(-log[10](p)),xlab="Position",type="n",cex.axis=1.5,cex.lab=1.5)
points(-log10(GEMMA_p_lrt)~Physical,data=chr10,col=cols,pch=20)
points(-log10(GEMMA_p_lrt)~Physical,data=caviar_snp,col="black",pch=1,cex=1.5,lwd=2)
axis(1,lwd=2,at=seq(1e5*floor(min(chr10$Physical)/1e5),1e5*ceiling(max(chr10$Physical)/1e5),0.1*1e6),labels=seq(0,700,100))
axis(2,las=1,lwd=2)
legend.col(col = colfunc(100), lev = chr10$R2)
dev.off()








