#! /usr/bin/Rscript --no-save --no-restore

library(LEA)

args=(commandArgs(TRUE))
#the first input is which environmental PC it is 
pc=args[1]
#the second input is how many latent group there is (K)
k=args[2]
print(pc) #PC1
print(k)  #1

LEA_summary<-function(pc=pc,k=k) {

wd=paste("/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/LEA/SwAsp.94samples_",pc,".lfmm/K",k,sep="")
setwd(wd)

z.table=NULL
for (i in 1:5) {
	file.name=paste(wd,"/run",i,"/SwAsp.94samples_r",i,"_s1.",k,".zscore",sep="")
	z.table=cbind(z.table,read.table(file.name)[,1])
}
z.score=apply(z.table,MARGIN=1,median)
lambda=median(z.score^2)/0.456
cat(lambda)
#lambda values:PC1: K1-0.450139, K2-0.4502006, K3-0.4510953

#ap.values=pchisq(z.score^2/lambda,df=1,lower=F)

#q=0.01
#L=length(ap.values)
#w=which(sort(ap.values)<q*(1:L)/L)
#candidates=order(ap.values)[w]

pos=read.table("/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/vcf/maf005/pos/all_snp.txt",header=T)

mean_pc_k=mean(z.score)
sd_pc_k=sd(z.score)

#number of SNPs
n=length(z.score)

p_value=NULL
for (i in 1:n)
{
if (z.score[i]>mean_pc_k)
{p_value[i]=pnorm(z.score[i],mean=mean_pc_k,sd=sd_pc_k,lower.tail=FALSE)}
else {p_value[i]=pnorm(z.score[i],mean=mean_pc_k,sd=sd_pc_k)}
}

summary_table=as.data.frame(cbind(pos$CHROM,pos$POS,z.score,p_value))
names(summary_table)=c("CHROM","POS","Z-score","P_value")

return(summary_table)
write.table(intercept_pc1,file=paste("LEA.",pc,".K",k,"z.score.txt"),sep="\t",quote=F,row.names=F,col.names=F)
}


