#! /usr/bin/Rscript --no-save --no-restore
.libPaths("/pica/h1/jingwang/R/x86_64-redhat-linux-gnu-library/3.2")

library(LDheatmap)
library(data.table)

args=(commandArgs(TRUE))
setwd(args[1])  ##working location
ld=read.table(args[2])  ##ld matrix
snp=read.table(args[3],header=F)  ##the example snp name

snp_name=snp$V1

ld_png=paste(args[2],".LDheatmap.png",sep="")
ld_pdf=paste(args[2],".LDheatmap.pdf",sep="")

#ld matrix is the r^2 matrix
ld=as.matrix(ld)
rownames(ld)=snp_name
colnames(ld)=snp_name

#colour of the ld plot 
#rgb.palette=colorRampPalette(rev(c("blue","yellow","red")),space="rgb")
rgb.palette=colorRampPalette(rev(c("grey90","red")),space="rgb")

#make the plot
#pdf(ld_pdf)
#png(ld_png,width = 6, height = 6, units = 'in', res=500)
pdf(file=ld_pdf,width = 6, height = 6)
#corrplot(cov,method="color",tl.col="black",tl.srt=45,addgrid.col="gray50",col=colorRampPalette(c("blue","white","red"))(100))
#LDheatmap(ld,add.map=FALSE,flip=T,color=rgb(20),title="")
LDheatmap(ld, SNP.name = c("Potra001246:25256"),color=rgb.palette(15),add.map=FALSE,title="")
#LDheatmap(ld, SNP.name = c("Potra001246:25256"),color=rgb.palette(15),add.map=TRUE,title="")
dev.off()


