#! /usr/bin/Rscript --no-save --no-restore
.libPaths("/pica/h1/jingwang/R/x86_64-redhat-linux-gnu-library/3.2")
library(data.table)

###the original output from GEMMA
gemma=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/budset/GEMMA/output/SwAsp.94.filter.all.budset.noPC.assoc.txt",header=T)
###calculate the Z-score of the GEMMA
gemma$z=gemma$beta/gemma$se

###the significant SNPs in chr10 region
chr10_sig=fread("/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/ld_r2_pos/sig_local.chr10.txt",header=T)

###make the gemma and chr10 as the same order in order to add the physical position to the gemma_chr10
gemma_chr10=gemma[which(gemma$rs %in% chr10_sig$SNP),]
gemma_chr10_order=gemma_chr10[order(gemma_chr10$rs),]
chr10_sig_order=chr10_sig[order(chr10_sig$SNP),]
gemma_chr10_order$chr10_pos=chr10_sig_order$chr10_pos
gemma_chr10_order_pos=gemma_chr10_order[order(gemma_chr10_order$chr10_pos),]

out=gemma_chr10_order_pos[,c("chr10_pos","z"),with=F]

write.table(out,file="chr10.caviar.txt",sep="\t", quote=F, row.names=F, col.names=F)


