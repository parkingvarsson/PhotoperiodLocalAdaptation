#! /bin/bash -l



plink="/proj/b2011141/tools/plink-1.07-x86_64/plink"
#plink="/proj/b2011141/tools/plink_linux_x86_64_dev/plink"
ped="/proj/b2011141/nobackup/all_populations/GATK/HC/snp_filter/GWAS/plink/asp201.SwAsp.gatk.hap.snp.rm_indel.bed.GQ10.biallelic.rm_miss_1.bed.plink.ped"
#map="/proj/b2011141/nobackup/all_populations/GATK/total/snp_filter/plink/asp201.chrom.map"
map="/proj/b2011141/nobackup/all_populations/GATK/HC/snp_filter/GWAS/plink/asp201.SwAsp.gatk.hap.snp.rm_indel.bed.GQ10.biallelic.rm_miss_1.bed.plink.map"


#$plink --ped $ped --map $map --indep-pairwise 50 5 0.2 
$plink --ped $ped --map $map --indep 50 5 2 
$plink --ped $ped --map $map --extract plink.prune.in --out pruneddata --recode






