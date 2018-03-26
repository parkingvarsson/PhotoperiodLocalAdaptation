#! /usr/bin/Rscript --no-save --no-restore
.libPaths("/pica/h1/jingwang/R/x86_64-redhat-linux-gnu-library/3.2")

#########################################################3
###Main purpose: the main aim of this script is to make both the manhanttan plot and qqplot for the output of the three analyses: PCAdapt, LFMM and GEMMA


library(data.table)
library("RColorBrewer")
source("/proj/b2011141/pipeline/R/local_adaptation_paper/manhattan/20171015/manhanttan_plot.R")

######set work directory
setwd("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/Overlap_three_methods")

#######read in the output data from the three files
summary=fread("gunzip -c PCAdapt.LEA.GEMMA.chr.summary.txt",header=T,stringsAsFactors=F)

#######remove four scaffolds within the selected region of chr10 that is likely to be wrongly mapped
scaff<-scan("/proj/b2011141/pipeline/R/local_adaptation_paper/manhattan/20171015/scaffold_likely_wrong.txt",what="char")
#[1] "Potra004011" "Potra000956" "Potra003230" "Potra158213"
grep(scaff[1],summary$SNP)->idx1
grep(scaff[2],summary$SNP)->idx2
grep(scaff[3],summary$SNP)->idx3
grep(scaff[4],summary$SNP)->idx4

summary<-summary[-c(idx1,idx2,idx3,idx4),]

#######read in the 1615 significnat SNPs we found in the end
sig_snp=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/Overlap_three_methods/snps/significant_snps.three.methods.txt",header=T)
scaffolds_7=c("Potra000799|Potra000908|Potra000342|Potra001246|Potra004002|Potra003230|Potra000530")
chr10_7scaffolds=summary[grepl(scaffolds_7,summary$SNP),]
sig_snp_chr10=chr10_7scaffolds$SNP

#######remove several another scaffolds which likely located in the chr10 region but are too short to contain in the 7 scaffolds
scaff2=c("Potra183460|Potra195127|Potra183319|Potra186780|Potra000530|Potra186900|Potra000554|Potra173873|Potra189005|Potra008516")
remove<-summary[grepl(scaff2,summary$SNP),]
summary=summary[-which(summary$SNP %in% remove$SNP),]

#######extract the output from the three results separately: PCAdapt, LFMM, GEMMA
##############################################################
#1.PCAdapt
pcadapt=summary[,c("CHROM","POS","SNP","PCAdapt_pvalue"),with=F]
names(pcadapt)=c("CHR","POS","SNP","P")
pcadapt_pvalue_cutoff=min(-log10(summary[which(summary$PCAdapt_qvalue<=0.05),PCAdapt_pvalue]))

#######make the plot
#pdf("PCAdapt.manhanttan.pdf",width=15,height=4,useDingbats=FALSE)
png("PCAdapt.manhanttan.png",width=15,height=4,units='in',res=500)
par(mar=c(5,5,2,2))
(pcadapt$P<0.01)->idx
#highlight=paste(sig_snp$CHROM,sig_snp$POS,sep=":")
highlight=sig_snp_chr10
manhattan(pcadapt[idx,],col =c("dodgerblue4","dodgerblue3"),highlight=highlight)
abline(h=pcadapt_pvalue_cutoff,lwd=2,lty=3,col="grey30")
dev.off()

##############################################################
#2. LFMM
lfmm=summary[,c("CHROM","POS","SNP","LEA_pvalue"),with=F]
names(lfmm)=c("CHR","POS","SNP","P")
lfmm_pvalue_cutoff=min(-log10(summary[which(summary$LEA_qvalue<=0.05),LEA_pvalue]))

#######make the plot
#pdf("LFMM.manhanttan.pdf",width=15,height=4,useDingbats=FALSE)
png("LFMM.manhanttan.png",width=15,height=4,units='in',res=500)
par(mar=c(5,5,2,2))
(lfmm$P<0.01)->idx
#highlight=paste(sig_snp_chr10$CHROM,sig_snp_chr10$POS,sep=":")
highlight=sig_snp_chr10
manhattan(lfmm[idx,],col =c("dodgerblue4","dodgerblue3"),highlight=highlight)
abline(h=lfmm_pvalue_cutoff,lwd=2,lty=3,col="grey30")
dev.off()

##############################################################
#3.GEMMA
gemma=summary[,c("CHROM","POS","SNP","GEMMA_pvalue"),with=F]
names(gemma)=c("CHR","POS","SNP","P")
gemma_pvalue_cutoff=min(-log10(summary[which(summary$GEMMA_qvalue<=0.05),GEMMA_pvalue]))

#######make the plot
#pdf("GEMMA.manhanttan.pdf",width=15,height=4,useDingbats=FALSE)
png("GEMMA.manhanttan.png",width=15,height=4,units='in',res=500)
par(mar=c(5,5,2,2))
(gemma$P<0.01)->idx
#highlight=paste(sig_snp_chr10$CHROM,sig_snp_chr10$POS,sep=":")
highlight=sig_snp_chr10
manhattan(gemma[idx,],col =c("dodgerblue4","dodgerblue3"),highlight=highlight)
abline(h=gemma_pvalue_cutoff,lwd=2,lty=3,col="grey30")
dev.off()


####################################################################################################
####################################################################################################
###2. Write the qq-plot function

ppoints<-function (n, a = if (n <= 10) 3/8 else 1/2)
{
    if (length(n) > 1L)
        n <- length(n)
    if (n > 0)
        (1L:n - a)/(n + 1 - 2 * a)
    else numeric()
}

qqplot<-function (pvector,critical,col, ...)
{
    if (!is.numeric(pvector))
        stop("Input must be numeric.")
    pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) &
        is.finite(pvector) & pvector < 1 & pvector > 0]
    o = -log10(sort(pvector, decreasing = FALSE))
    e = -log10(ppoints(length(pvector)))
    c = -log10(critical)
    o_e=data.frame(cbind(o,e))
    o_outlier=o_e[which(o_e$o>c),]
    o_expected=o_e[which(o_e$o<=c),]
    plot(x=o_expected$e,y=o_expected$o,pch=20,cex.lab=2,cex.axis=1.5,axes=F,xlim = c(0, max(e)), ylim = c(0,max(o)), xlab = expression(Expected ~ ~-log[10](italic(p))),ylab = expression(Observed ~ ~-log[10](italic(p))))
    axis(1)
    axis(2)
    par(new=T)
    points(o_outlier$e,o_outlier$o,pch=20,col=col)
    abline(0, 1, col = "black")
#    abline(h=c,col=col,lty=3,lwd=2)
}

########read in the data
data=fread("PCAdapt.LEA.GEMMA.summary.txt",header=T)
######################################################################################################
###1 Remove SNPs showing the excess of heterzygotes, HWE<1e-8
###The SNPs with HWE-pvalue <1e-8 is located in /proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/maf005/pos/hardy/hwe.extreme.pos
hwe_pos=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/maf005/pos/hardy/hwe.extreme.pos",header=F)

data$SNP=paste(data$CHROM,data$POS,sep=":")
hwe_pos$SNP=paste(hwe_pos$V1,hwe_pos$V2,sep=":")
data_new=data[which(!data$SNP %in% hwe_pos$SNP),]


##2.1pcadapt
png(filename="pcadapt.qqplot.png",width=4,height=4,units='in',res=500)
#pdf("pcadapt.qqplot.pdf",width=5,height=5)
par(mfrow=c(1,1))
par(mar=c(5,5,1,1))
##qqplot for PCAdapt
pcadapt_critical=max(data_new[which(data_new$PCAdapt_qvalue<0.05),PCAdapt_pvalue])
qqplot(data_new$PCAdapt_pvalue,pcadapt_critical,col="firebrick1")
dev.off()

##2.2 LFMM
png(filename="lea.qqplot.png",width=4,height=4,units='in',res=500)
#pdf("lea.qqplot.pdf",width=5,height=5)
par(mfrow=c(1,1))
par(mar=c(5,5,1,1))
LEA_critical=max(data_new[which(data_new$LEA_qvalue<0.05),LEA_pvalue])
qqplot(data_new$LEA_pvalue,LEA_critical,col="firebrick1")
dev.off()

##2.3 GEMMA
png(filename="gemma.qqplot.png",width=4,height=4,units='in',res=500)
#pdf("gemma.qqplot.pdf",width=5,height=5)
par(mfrow=c(1,1))
par(mar=c(5,5,1,1))
GEMMA_critical=max(data[which(data$GEMMA_qvalue<0.05),GEMMA_p_lrt])
qqplot(data_new$GEMMA_p_lrt,GEMMA_critical,col="firebrick1")
dev.off()



