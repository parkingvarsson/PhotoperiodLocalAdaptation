#! /usr/bin/Rscript --no-save --no-restore
.libPaths("/pica/h1/jingwang/R/x86_64-redhat-linux-gnu-library/3.2")

library(data.table)
library("RColorBrewer")
source("/proj/sllstore2017050/pipeline_jingwang/pipeline/R/local_adaptation_paper/manhattan/20171015/legend.col.R")

####################################################################
###Main aim: to create the GEMMA plot specifically for the Chr10 selected region, and also include the LD pattern (with the most significant SNP: Potra001245:25256) in it
setwd("/proj/sllstore2017050/nobackup/milou_files/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/sliding_three_methods")

chr10=fread("chr10.7scaffolds.region.three.txt",header=T)
ld<-fread("/proj/sllstore2017050/nobackup/milou_files/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/ld_r2_pos/SwAsp.chr10.7scaffolds.recode.Potri.rm_het.gt.beagle.hap.Dprimer.ld",header=T)

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

##subset to the last 30kb of Potra001246
chr10[grep("Potra001246",chr10$SNP),]->chr10
chr10[chr10$Physical>377589, ]->chr10
ft2<-fread("/proj/sllstore2017050/nobackup/milou_files/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/sliding_three_methods/FT2_gene/FT2_transcript.txt",head=T)

pdf("GWAS_FT2_Dprime.pdf",width=5,height=6)
par(mgp=c(2.25,1,0))
colfunc <- colorRampPalette(c("dodgerblue4","darkorange","firebrick1"))
cols<-colfunc(10)[ceiling(abs(chr10$Dprime)*10)]
plot(-log10(GEMMA_p_lrt)~ Physical,data=chr10,yaxt="n",xaxt="n",bty="n",ylab=expression(-log[10](p)),xlab="Position",type="n",cex.axis=1.5,cex.lab=2,ylim=c(0,28.5))
points(-log10(GEMMA_p_lrt)~Physical,data=chr10,col=cols,pch=20)
points(-log10(GEMMA_p_lrt)~Physical,data=caviar_snp,col="black",pch=1,cex=1.5,lwd=2)
axis(1,lwd=2,at=c(380000,390000,400000,410000,420000),labels=c(380,390,400,410,420),cex.axis=1.5)
axis(2,las=1,lwd=2,cex.axis=1.5)
for(i in 1:dim(ft2)[1]){
	if(ft2[i,"gene",with=F]=="FT2"|ft2[i,"gene",with=F]=="FT2beta"){
		col<-ifelse(i %in% grep("utr",ft2$annotation),"light blue","red")
		rect(xleft=421484-ft2[i,"stop",with=F],xright=421484-ft2[i,"start",with=F],ybottom=25,ytop=26,col=col,border=NA)
	}else{
		col<-ifelse(i %in% grep("utr",ft2$annotation),"grey75","grey15")
		rect(xleft=421484-ft2[i,"stop",with=F],xright=421484-ft2[i,"start",with=F],ybottom=25,ytop=26,col=col,border=NA)	
	}
}
arrows(x1=421484-max(ft2[gene=="FT2","stop",with=F]),x0=421484-min(ft2[gene=="FT2","start",with=F]),y0=26.5,y1=26.5,lwd=2,length=0.1)
arrows(x1=421484-max(ft2[gene=="FT2beta","stop",with=F]),x0=421484-min(ft2[gene=="FT2beta","start",with=F]),y0=26.5,y1=26.5,lwd=2,length=0.1)
arrows(x1=421484-max(ft2[gene=="INO80","stop",with=F]),x0=421484-min(ft2[gene=="INO80","start",with=F]),y0=26.5,y1=26.5,lwd=2,length=0.1)
arrows(x1=421484-max(ft2[gene=="FAS1","stop",with=F]),x0=421484-min(ft2[gene=="FAS1","start",with=F]),y0=26.5,y1=26.5,lwd=2,length=0.1)
arrows(x1=421484-max(ft2[gene=="AAD5","stop",with=F]),x0=421484-min(ft2[gene=="AAD5","start",with=F]),y0=26.5,y1=26.5,lwd=2,length=0.1)
#text(x=382000,y=27.5,"PtAAD5 (Potra001246g10696)",cex=1)
text(x=380000,y=27.5,"PtAAD5",cex=1.2)
#text(x=391000,y=29,"PtFAS1 (Potra001246g10695)",cex=1)
text(x=391000,y=29,"PtFAS1",cex=1.2)
#text(x=397000,y=27.5,"PtFT2 (Potra001246g10694)",cex=1)
text(x=397000,y=27.5,"PtFT2",cex=1.2)
#text(x=409000,y=29,"PtINO80 (Potra001246g10693)",cex=1)
text(x=409000,y=29,"PtINO80",cex=1.2)
text(x=418000,y=27.5,expression(paste("PtFT2",beta)),cex=1.2)
legend.col(col = colfunc(100), lev = chr10$Dprime)
dev.off()



