## Create imputing used .ped file with different level of missing value ### 
rm(list=ls())
setwd("~/Documents/populus")
## another way to do the above procedure which has been done by Perl

simulate <- read.table("miss_first100s_re.ped", sep=" ",header=F)
title_ped <- simulate[,1:6]
simulate_m <- simulate[,7:ncol(simulate)]
simulate_m <- as.matrix(simulate_m)
total_gene <- dim(simulate_m)[1]*dim(simulate_m)[2]

missing <- c(0.05, 0.1,0.15,0.2,0.25,0.3,0.4,0.5) # setting a missing percent

k=0.05
for (k in missing){
    #set.seed(1234)
    missNum <- round(total_gene*k)   # the # to be replaced as missing value
    miss <- sample(total_gene,missNum,replace=F) # random select the missing genotype
    simu_miss <- simulate_m  # creat a copy of genotype
    simu_miss[miss] <- "0/0" # replace the genotype, data should be as a matrix
    re_miss <- gsub(pattern="/", replacement=" ",simu_miss)
    simu <- data.frame(cbind(title_ped, re_miss), row.names=NULL)
    name <- paste("simu_", k, "_miss_first100s_try.ped", sep="")
    write.table(simu,file=name,sep=" ",row.names=F,col.names=F, quote=F)
}

