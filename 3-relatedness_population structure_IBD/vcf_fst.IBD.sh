#! /bin/bash -l



module load bioinfo-tools
#module load vcftools

vcftools="/proj/b2011141/tools/vcftools_0.1.12/bin/vcftools"
bgzip="/proj/b2011141/tools/bgzip"

Inputvcf=$1
VCFDir=`dirname $1`

Out=${Inputvcf##*/}
echo $Out

fst="/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/IBD"

OutDir=$fst/$2_$3

if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

pop="/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/fst/pop"

pop1=$pop/$2.txt
pop2=$pop/$3.txt

Out_fst=${Out%.vcf.gz}.$2.$3

$vcftools --gzvcf $1 --weir-fst-pop $pop1 --weir-fst-pop $pop2 --out $OutDir/$Out_fst



