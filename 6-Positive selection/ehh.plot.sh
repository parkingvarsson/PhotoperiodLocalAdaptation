#! /bin/bash -l


module load bioinfo-tools
#module load R
module load selscan

colormap="/proj/b2011141/pipeline/R/selscan/colormap.plotting.R"
tped="/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/selscan/chr10/chr10.7scaffolds.ihs.tped"
OutDir="/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/selscan/chr10/ehh/"

snp=$1
##split snp name into scaffold name and snp position
#read scaffold <<< $(echo $1 | awk '{split($0,a,":"); print a[1]}')
#echo $scaffold

#selscan --ehh $snp --tped ${scaffold}.tped --ehh-win 50000
#selscan --ehh $snp --tped $tped --out $OutDir/outfile --ehh-win 50000
#selscan --ehh $snp --tped $tped --out $OutDir/outfile --ehh-win 15000
Rscript $colormap $OutDir $snp


