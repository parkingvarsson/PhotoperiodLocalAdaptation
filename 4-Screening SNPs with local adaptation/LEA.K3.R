#! /usr/bin/Rscript --no-save --no-restore

setwd("/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/LEA")

library(LEA)

args=(commandArgs(TRUE))
print(args)
env=paste(args,".env",sep="")
project=lfmm("SwAsp.94samples.lfmm",env,K=3,CPU=8,repetitions=5,project="continue")


