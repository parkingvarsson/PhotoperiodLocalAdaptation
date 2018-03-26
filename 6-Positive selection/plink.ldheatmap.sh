#! /bin/bash -l


module load bioinfo-tools
module load plink
###Main aim: simply script to use norm function from selscan 

ld_test=$1  ###choose the step
group=$2  ###north, mid, south, all

InputDir="/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/plink"
tped_r="/proj/b2011141/pipeline/R/local_adaptation_paper/LDheatmap/chr10.sig_snps.tped.R"  ##extracting significant SNPs in chr10
tped_ft2_r="/proj/b2011141/pipeline/R/local_adaptation_paper/LDheatmap/ft2_gene.tped.R" ##extracting SNPs around FT2 gene (up and downstream 2kb)
ld_dir=$InputDir/ld

if [ ! -d "$ld_dir" ]; then
mkdir -p $ld_dir
fi

tped=$InputDir/SwAsp.chr10.7scaffolds.$group.tped
tfam=$InputDir/SwAsp.chr10.7scaffolds.$group.tfam

########the first step is to calculate the LD r2 matrix correlation
if [ "$ld_test" == "1" ]; then

#plink --tped $tped --tfam $tfam --r2 square dprime --allow-no-sex --out $ld_dir/SwAsp.chr10.7scaffolds.$group --distance-matrix
plink --tped $tped --tfam $tfam --r2 dprime --allow-no-sex --out $ld_dir/SwAsp.chr10.7scaffolds.$group --distance-matrix

elif [ "$ld_test" == "2" ]; then  ###use Rscript to extract significant SNPs and then use plink to create pairwise LD matrix

Rscript $tped_r $InputDir SwAsp.chr10.7scaffolds.$group.tped
plink --tped $tped.sig.tped --tfam $tfam --r2 square --allow-no-sex --out $ld_dir/SwAsp.chr10.7scaffolds.$group.sig_snps --distance-matrix

elif [ "$ld_test" == "3" ]; then  ###use Rscript to extract SNPs around FT2 gene and then use plink to create pairwise LD matrix

Rscript $tped_ft2_r $InputDir SwAsp.chr10.7scaffolds.$group.tped
plink --tped $tped.ft2_2kb.tped --tfam $tfam --r2 square --allow-no-sex --out $ld_dir/SwAsp.chr10.7scaffolds.$group.ft2 --distance-matrix
plink --tped $tped.ft2_2kb.sig_snps.tped --tfam $tfam --r2 square --allow-no-sex --out $ld_dir/SwAsp.chr10.7scaffolds.$group.ft2.sig_snps --distance-matrix

elif [ "$ld_test" == "4" ]; then

ld_plot="/proj/b2011141/pipeline/R/local_adaptation_paper/LDheatmap/LDheatmap.R"
Rscript $ld_plot $ld_dir $ld_dir/SwAsp.chr10.7scaffolds.$group.ld $ld_dir/snp.name

elif [ "$ld_test" == "5" ]; then

ld_plot="/proj/b2011141/pipeline/R/local_adaptation_paper/LDheatmap/LDheatmap.R"
Rscript $ld_plot $ld_dir $ld_dir/SwAsp.chr10.7scaffolds.$group.sig_snps.ld $ld_dir/sig.snp.name

elif [ "$ld_test" == "6" ]; then

ld_plot="/proj/b2011141/pipeline/R/local_adaptation_paper/LDheatmap/LDheatmap.R"
Rscript $ld_plot $ld_dir $ld_dir/SwAsp.chr10.7scaffolds.$group.ft2.ld $ld_dir/ft2.snp.name
Rscript $ld_plot $ld_dir $ld_dir/SwAsp.chr10.7scaffolds.$group.ft2.sig_snps.ld $ld_dir/ft2.sig.snp.name

fi
