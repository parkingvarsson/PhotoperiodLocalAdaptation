#! /usr/bin/Rscript --no-save --no-restore


###This script is used to estimate the genetic distances between populations

setwd("/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/IBD/")


###step1, calculate the pairwise fst between populations
n=12
cat(n)
d_genetic=matrix(NA,n,n)
for (i in 1:11) {
	b=i+1
	for (j in b:12) {
		file=paste("pop",i,"_pop",j,sep="")	
		filename=paste("/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/IBD/",file,"/SwAsp_94samples.filter.maf005.pop",i,".pop",j,".weir.fst",sep="")
		cat(filename)
		fst=read.table(filename,header=T)
		fst_mean=mean(fst$WEIR_AND_COCKERHAM_FST,na.rm=T)
		cat(fst_mean)
		d_genetic[i,j]=fst_mean
		}
	}
###step2, for the matrix table, we also need to fill the other part of the table, and also fill the diagonal by 0
for (i in 1:11) {
	c=i+1
	for (j in c:12) {
		d_genetic[j,i]=d_genetic[i,j]
	}
}
for (i in 1:12) {
	d_genetic[i,i]=0
}

write.table(d_genetic,file="genetic-distances.txt",sep="\t",col.names=FALSE, row.names=FALSE, quote=FALSE)


