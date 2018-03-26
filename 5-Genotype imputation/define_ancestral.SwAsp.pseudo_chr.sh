#! /bin/bash -l

#SBATCH -A b2010014
#SBATCH -p core
#SBATCH -n 1
#SBATCH -o define_ancestral.pseudo.out
#SBATCH -e define_ancestral.pseudo.err
#SBATCH -J define_ancestral.pseudo.job
#SBATCH -t 12:00:00

module load bioinfo-tools
module load vcftools
module load samtools

###Step0: Before all downstream steps, I think the SNPs with high heterzygosity (and with extreme P-values from Hardy-Weinberg test) should be removed since they confould the haplotype homozygosity singals, especially in FT2 region


###Step1: extract the SNPs of chromosome 10 from the Inputvcf file 
Inputvcf="/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/maf005/pos/SwAsp_94samples.filter.maf005.pos.recode.Potri.vcf.gz"

vcf_Dir=`dirname $Inputvcf`

####################################################################################################
#step1.1  Before the step3, I need to remove the SNPs showing extreme P_HET_EXCESS from HWE 

rm_het_hwe=$vcf_Dir/rm_het_hwe

if [ ! -d "$rm_het_hwe" ]; then
mkdir -p $rm_het_hwe
fi

hwe="$vcf_Dir/hardy/SwAsp_94samples.filter.maf005.pos.recode.Potri.hwe"

#awk '$8<0.00000001' $hwe | cut -f 1,2 > $rm_het_hwe/all.hwe.rm.pos
###remove those SNPs from the pseudo-chromsome datasets
#vcftools --gzvcf $Inputvcf --exclude-positions $rm_het_hwe/all.hwe.rm.pos --recode --recode-INFO-all --out $rm_het_hwe/SwAsp_94samples.filter.maf005.pos.Potri.rm_het
#bgzip $rm_het_hwe/SwAsp_94samples.filter.maf005.pos.Potri.rm_het.recode.vcf
###remove those SNPs from original scaffold dataset
#vcftools --gzvcf $vcf_Dir/SwAsp_94samples.filter.maf005.pos.recode.vcf.gz --exclude-positions $rm_het_hwe/all.hwe.rm.pos --recode --recode-INFO-all --out $rm_het_hwe/SwAsp_94samples.filter.maf005.pos.rm_het
#bgzip $rm_het_hwe/SwAsp_94samples.filter.maf005.pos.rm_het.recode.vcf


###Step2:use BEAGLE to do imputation
beagle="/proj/b2011141/tools/beagle4.1/beagle.27Jul16.86a.jar"

beagle_dir=$rm_het_hwe/beagle

if [ ! -d "$beagle_dir" ]; then
mkdir -p $beagle_dir
fi

module load java

Out="SwAsp_94samples.filter.maf005.pos.Potri.rm_het.recode.vcf.gz"
Out_beagle=${Out%.recode.vcf.gz}.gt
Out_beagle_phase=${Out%.recode.vcf.gz}.gt.beagle
##the first is to impute the missing genotype
#java -Xmx100g -jar $beagle gtgl=$rm_het_hwe/$Out out=$beagle_dir/$Out_beagle nthreads=12
#java -Xmx100g -jar $beagle gt=$beagle_dir/${Out_beagle}.vcf.gz out=$beagle_dir/$Out_beagle_phase nthreads=12

##transfer to plink file
vcf_plink="/proj/b2011141/pipeline/vcftools/SwAsp/vcf_plink.sh"

$vcf_plink $beagle_dir/${Out_beagle_phase}.vcf.gz
zcat $beagle_dir/${Out_beagle_phase}.vcf.gz |grep -v "#" |cut -f 1 |sed 's/PseudoChr0//g' |sed 's/PseudoChr//g' > $beagle_dir/plink/chr
cut -f 2- $beagle_dir/plink/${Out_beagle_phase}.tped > $beagle_dir/plink/temp 
paste $beagle_dir/plink/chr $beagle_dir/plink/temp > $beagle_dir/plink/${Out_beagle_phase}.tped && rm $beagle_dir/plink/chr $beagle_dir/plink/temp

##-------Step4: based on the information of ancestral state for SNPs, transfer the tped to ihs.tped where 0 represents ancestral allele and 1 represents derived allele

tped_anc="/proj/b2011141/pipeline/perl/SwAsp/tped_anc.01.pl"
anc_Dir="/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/maf005/pos/ancestral"
perl $tped_anc $beagle_dir/plink/${Out_beagle_phase}.tped  $anc_Dir/SwAsp.pos.outgroup.final.anc




