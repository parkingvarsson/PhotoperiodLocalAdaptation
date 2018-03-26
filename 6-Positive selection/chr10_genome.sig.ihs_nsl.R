#! /usr/bin/Rscript --no-save --no-restore
.libPaths("/pica/h1/jingwang/R/x86_64-redhat-linux-gnu-library/3.2")
library(RColorBrewer)
library(dplyr)
colors <- brewer.pal(10,"Paired")

library(data.table)

setwd("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/selscan/sig_test")

#####To test for whether there is significant concentration of selection signals on the region surrounding the PtFT2, we divided the 19 pseudo-chromosomes (without the seven scaffolds around the PtFT2 locus) into non-overlapping windows of 700 kb and calculated the proportion of SNPs with |iHS| > 2 or with [nSL|>2 in each window. 


####function
#ihs_mean
slideMean_ihs<-function(x,chr,windowsize=700000,slide=700000){
  end_pos=x[nrow(x),]$POS
  idx1<-seq(1,end_pos,by=slide);
  idx1+windowsize->idx2;
  new_pos=(idx2-idx1)/2+idx1
  
  ihs_2_n=c()
  snp_n=c()
  ihs_2_freq=c()
  chrom=c()
  for (i in 1:length(new_pos)){
    chrom[i]=chr
    ihs_2_n[i]=length(which(x[which(x$POS<idx2[i] & x$POS>idx1[i]),]$n_ihs_2==1))
    snp_n[i]=nrow(x[which(x$POS<idx2[i] & x$POS>idx1[i]),])
    ihs_2_freq[i]=length(which(x[which(x$POS<idx2[i] & x$POS>idx1[i]),]$n_ihs_2==1))/nrow(x[which(x$POS<idx2[i] & x$POS>idx1[i]),])
    }
  
  ihs_window=data.frame(cbind(chrom,new_pos,ihs_2_n,snp_n,ihs_2_freq))
  return(ihs_window)
}


genome_ihs=data.frame()

for (chr in 1:19){
ihs_chr=paste("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/maf005/pos/rm_het_hwe/selscan/ihs/SwAsp.genome.all.chr",chr,".ihs.out.100bins.norm",sep="")
ihs=fread(ihs_chr,header=F)
names(ihs)=c("SNP","POS","Derived_Freq","iHH1","iHH0","iHS","norm_iHS","n_ihs_2")
ihs_w=slideMean_ihs(ihs,chr)
genome_ihs=rbind(genome_ihs,ihs_w)
}

chr10_ihs=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/selscan/ihs/SwAsp.chr10.7scaffolds.all.ihs.out.100bins.norm",header=F)
names(chr10_ihs)=c("SNP","POS","Derived_Freq","iHH1","iHH0","iHS","norm_iHS","n_ihs_2")
chr10_ihs_2_freq=length(which(chr10_ihs$n_ihs_2=="1"))/nrow(chr10_ihs)

pdf(file="chr10.ihs_sig.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
hist(genome_ihs[which(genome_ihs$snp_n>1000),]$ihs_2_freq,breaks=50,xlim=c(0,0.23),col="grey50",xlab="Prop of SNPs with |iHS|>2",main="",cex.lab=1.5)
abline(v=chr10_ihs_2_freq,col="dark red",lwd=2)  ###ranking P value is 1/627:0.001594896
dev.off()

##nsl_function
slideMean_nsl<-function(x,chr,windowsize=700000,slide=700000){
  end_pos=x[nrow(x),]$POS
  idx1<-seq(1,end_pos,by=slide);
  idx1+windowsize->idx2;
  new_pos=(idx2-idx1)/2+idx1

  nsl_2_n=c()
  snp_n=c()
  nsl_2_freq=c()
  chrom=c()
  for (i in 1:length(new_pos)){
    chrom[i]=chr
    nsl_2_n[i]=length(which(x[which(x$POS<idx2[i] & x$POS>idx1[i]),]$n_nsl_2==1))
    snp_n[i]=nrow(x[which(x$POS<idx2[i] & x$POS>idx1[i]),])
    nsl_2_freq[i]=length(which(x[which(x$POS<idx2[i] & x$POS>idx1[i]),]$n_nsl_2==1))/nrow(x[which(x$POS<idx2[i] & x$POS>idx1[i]),])
    }

  nsl_window=data.frame(cbind(chrom,new_pos,nsl_2_n,snp_n,nsl_2_freq))
  return(nsl_window)
}


genome_nsl=data.frame()

for (chr in 1:19){
nsl_chr=paste("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/maf005/pos/rm_het_hwe/selscan/nsl/SwAsp.genome.all.chr",chr,".nsl.out.100bins.norm",sep="")
nsl=fread(nsl_chr,header=F)
names(nsl)=c("SNP","POS","Derived_Freq","sL1","sL0","nSL","norm_nSL","n_nsl_2")
nsl_w=slideMean_nsl(nsl,chr)
genome_nsl=rbind(genome_nsl,nsl_w)
}

chr10_nsl=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/selscan/nsl/SwAsp.chr10.7scaffolds.all.nsl.out.100bins.norm",header=F)
names(chr10_nsl)=c("SNP","POS","Derived_Freq","sL1","sL0","nSL","norm_nSL","n_nsl_2")
chr10_nsl_2_freq=length(which(chr10_nsl$n_nsl_2=="1"))/nrow(chr10_nsl)

pdf(file="chr10.nsl_sig.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
hist(genome_nsl[which(genome_nsl$snp_n>1000),]$nsl_2_freq,breaks=50,xlim=c(0,0.2),col="grey50",xlab="Prop of SNPs with |nSL|>2",main="",cex.lab=1.5)
##p-value=length(which(genome_nsl[which(genome_nsl$snp_n>1000),]$nsl_2_freq>chr10_nsl_2_freq))/nrow(genome_nsl)
abline(v=chr10_nsl_2_freq,col="dark red",lwd=2)  ###ranking P value is 23/627=0.036682
dev.off()











