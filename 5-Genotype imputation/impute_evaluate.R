# This script is for evaluating the imputation accuracy in populus dataset
rm(list=ls())
setwd("~/Documents/populus/Aspen/imputation/beagle4.1/chr2_iter100")
MAF <- read.table("../noTitle_ch2.frq",header=T,stringsAsFactors = F)
miss_5 <- read.table("0.05_SNP_error_rate_ADD.txt",header=T,stringsAsFactors=F)
miss_5 <- miss_5[1:nrow(miss_5)-1,]
miss_10 <- read.table("0.1_SNP_error_rate_ADD.txt",header=T,stringsAsFactors=F)
miss_10 <- miss_10[1:nrow(miss_10)-1,]
miss_15 <- read.table("0.15_SNP_error_rate_ADD.txt",header=T,stringsAsFactors=F)
miss_15 <- miss_15[1:nrow(miss_15)-1,]
miss_20 <- read.table("0.2_SNP_error_rate_ADD.txt",header=T,stringsAsFactors=F)
miss_20 <- miss_20[1:nrow(miss_20)-1,]
miss_25 <- read.table("0.25_SNP_error_rate_ADD.txt",header=T,stringsAsFactors=F)
miss_25 <- miss_25[1:nrow(miss_25)-1,]
miss_30 <- read.table("0.3_SNP_error_rate_ADD.txt",header=T,stringsAsFactors=F)
miss_30 <- miss_30[1:nrow(miss_30)-1,]
miss_40 <- read.table("0.4_SNP_error_rate_ADD.txt",header=T,stringsAsFactors=F)
miss_40 <- miss_40[1:nrow(miss_40)-1,]
miss_50 <- read.table("0.5_SNP_error_rate_ADD.txt",header=T,stringsAsFactors=F)
miss_50 <- miss_50[1:nrow(miss_50)-1,]
accuracy <- data.frame(SNP=MAF$SNP,MAF=MAF$MAF,accu_05=1-miss_5$ratio,accu_10=1-miss_10$ratio,
                       accu_15=1-miss_15$ratio,accu_20=1-miss_20$ratio,accu_25=1-miss_25$ratio,
                       accu_30=1-miss_30$ratio,accu_40=1-miss_40$ratio,accu_50=1-miss_50$ratio,stringsAsFactors = F)
sample <- sample(x = 1:nrow(accuracy),size=100000,replace=F)
accuracy <- accuracy[sample,]
## 1. create MAF vs. accuracy in different missing level ##

#accuracy <- false
#accuracy[,6:ncol(accuracy)] <- 1-accuracy[,6:ncol(accuracy)] #change error rate to correct rate
head(accuracy)
## Fiting a curve line
acc_order <- accuracy[order(accuracy$MAF),]
rownames(acc_order)<- NULL
beagle_5 <- loess(accu_05~MAF, data=acc_order)
beagle_10 <- loess(accu_10~MAF, data=acc_order)
beagle_15 <- loess(accu_15~MAF, data=acc_order)
beagle_20 <- loess(accu_20~MAF, data=acc_order)
beagle_25 <- loess(accu_25~MAF, data=acc_order)
beagle_30 <- loess(accu_30~MAF, data=acc_order)
beagle_40 <- loess(accu_40~MAF, data=acc_order)
beagle_50 <- loess(accu_50~MAF, data=acc_order)
#beagle_50 <- loess(r_50_beagle~MAF, data=acc_order)

shapeit_5 <- loess(r_5_shapeit~MAF, data=acc_order)
shapeit_10 <- loess(r_10_shapeit~MAF, data=acc_order)
shapeit_15 <- loess(r_15_shapeit~MAF, data=acc_order)
shapeit_20 <- loess(r_20_shapeit~MAF, data=acc_orsder)
shapeit_25 <- loess(r_25_shapeit~MAF, data=acc_order)
shapeit_30 <- loess(r_30_shapeit~MAF, data=acc_order)
shapeit_35 <- loess(r_35_shapeit~MAF, data=acc_order)
shapeit_40 <- loess(r_40_shapeit~MAF, data=acc_order)
shapeit_50 <- loess(r_50_shapeit~MAF, data=acc_order)
# Draw a curve line indicate the MAF vs. accuracy
pdf(file = "AccuracyVSMAF.pdf",width = 6,height = 6)
par(mar=c(4,4,1.5,1.5))
plot(accu_15~MAF, data=acc_order, type="n",ylim=c(0.94,1), ylab="Imputation Accuracy",
     xlab="Minor Allele Frequency") #main="MAF vs. Imputation Accuracy", 
lines(acc_order$MAF, predict(beagle_5), col = "brown1", lwd=3)
lines(acc_order$MAF, predict(beagle_10), col = "deepskyblue", lwd=3)
lines(acc_order$MAF, predict(beagle_15), col = "darkorchid1", lwd=3)
lines(acc_order$MAF, predict(beagle_20), col = "dodgerblue3", lwd=3)
lines(acc_order$MAF, predict(beagle_25), col = "chartreuse4", lwd=3)
lines(acc_order$MAF, predict(beagle_30), col = "palevioletred2", lwd=3)
#lines(acc_order$MAF, predict(beagle_35), col = "darkgreen", lwd=2)
lines(acc_order$MAF, predict(beagle_40), col = "orange", lwd=3)
lines(acc_order$MAF, predict(beagle_50), col = "royalblue4", lwd=3)
legend(0.24,1, c("5% missing","10% missing", "15% missing","20% missing"),cex=0.9,
       lty=c(1,1,1,1,1),col=c("brown1","deepskyblue","darkorchid1","dodgerblue3"),
       bty="n",text.width =0.5,lwd=2)
legend(0.37,1, c("25% missing","30% missing","40% missing","50% missing"),cex=0.9,
       lty=c(1,1,1,1,1),col=c("chartreuse4","palevioletred2","orange","royalblue4"),
       bty="n",text.width =0.5,lwd=2)
dev.off()

par(mar=c(4,4,2.5,1.5))
plot(r_5_shapeit~MAF, data=acc_order, type="n",ylim=c(0.6,1), ylab="Imputation Accuracy",
     main="MAF vs. Imputation Accuracy by Shapeit", xlab="Minor Allele Frequency")
lines(acc_order$MAF, predict(shapeit_5), col = "brown1", lwd=2)
lines(acc_order$MAF, predict(shapeit_10), col = "deepskyblue", lwd=2)
lines(acc_order$MAF, predict(shapeit_15), col = "darkorchid1", lwd=2)
lines(acc_order$MAF, predict(shapeit_20), col = "dodgerblue3", lwd=2)
lines(acc_order$MAF, predict(shapeit_25), col = "chartreuse4", lwd=2)
lines(acc_order$MAF, predict(shapeit_30), col = "palevioletred2", lwd=2)
lines(acc_order$MAF, predict(shapeit_35), col = "darkgreen", lwd=2)
lines(acc_order$MAF, predict(shapeit_40), col = "orange", lwd=2)
lines(acc_order$MAF, predict(shapeit_50), col = "royalblue4", lwd=2)
legend(0.28,0.7, c("5% miss","10% miss", "15% miss","20% miss","25% miss"),cex=0.7,
       lty=c(1,1,1,1,1),col=c("brown1","deepskyblue","darkorchid1","dodgerblue3","chartreuse4"),
       bty="n",text.width =0.5,lwd=2)
legend(0.39,0.7, c("30% miss","35% miss", "40% miss","45% miss"),cex=0.7,
       lty=c(1,1,1,1,1),col=c("palevioletred2","darkgreen","orange","royalblue4"),
       bty="n",text.width =0.5,lwd=2)

