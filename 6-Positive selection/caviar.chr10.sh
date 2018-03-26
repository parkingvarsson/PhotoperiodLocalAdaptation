#! /bin/bash -l



caviar="/proj/b2011141/tools/caviar-master/CAVIAR-C++/CAVIAR"
caviar_input_r="/proj/b2011141/pipeline/local_adaptation_paper/caviar/caviar.z-score.R"

InputDir="/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/caviar"

###r script to create the chr10.caviar.txt
#Rscript $caviar_input_r
###use caviar to choose the the "causal" SNPs
##use R2 as LD measurement
#$caviar -o $InputDir/caviar.output.txt -z $InputDir/chr10.caviar.txt -l $InputDir/chr10.LD.r2.txt -r 0.99
##use Dprime as LD measurement
$caviar -o $InputDir/caviar.output.Dprime.txt -z $InputDir/chr10.caviar.txt -l $InputDir/chr10.sig_snps.LD.D.txt -r 0.99


