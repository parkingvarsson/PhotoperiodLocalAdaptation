#! /bin/bash -l


module load bioinfo-tools
module load python/2.7.11
###Main aim: the main aim here is to estimate a number of selection tests on the region around chr10: diversity, Fst, dxy, iHS, nSL, H12, SweepFinder-CLR
Inputvcf="/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/SwAsp.chr10.7scaffolds.recode.Potri.rm_het.gt.beagle.vcf.gz"

Out=${Inputvcf##*/}
echo $Out

vcf_Dir=`dirname $Inputvcf`

###################################################################################################################
##I skip this step now, since I found the results looks weird to me
#####Step1: estimate pi, Fst and dxy using the python script

###Step1.1:Parsing VCF files, use python script to transfer vcf file into .geno file, with both un-phased data and phased data
parse_vcf="/proj/b2011141/tools/twisst/genomics_general-master/VCF_processing/parseVCF.py"
###Step2.1 parsing un-phhased .vcf dataset
#summary_dir=$vcf_Dir/diversity_fst_dxy

#if [ ! -d "$summary_dir" ]; then
#mkdir -p $summary_dir
#fi

#python $parse_vcf -i $Inputvcf |gzip > $summary_dir/${Out%.vcf.gz}.geno.gz


###Step1.2: Diversity and divergence analyses in sliding windows separately for individuals from northern, middle and southern populations
popgen_window="/proj/b2011141/tools/twisst/genomics_general-master/popgenWindows.py"

###extract individual names in each groups of populations
#south_ind=$(tr '\n' ',' < $vcf_Dir/south.inds | sed 's/\,$//g')
#mid_ind=$(tr '\n' ',' < $vcf_Dir/mid.inds | sed 's/\,$//g')
#north_ind=$(tr '\n' ',' < $vcf_Dir/north.inds | sed 's/\,$//g')

#echo $south_ind
#echo $mid_ind
#echo $north_ind

###define widnow size: first try 5kb with 2kb sliding windows

#python $popgen_window -w 5000 -s 2500 -o 2500 -m 10 --windType coordinate -g $summary_dir/${Out%.vcf.gz}.geno.gz -o $summary_dir/${Out%.vcf.gz}.5kb_2kb.popgen.csv.gz -p S $south_ind -p M $mid_ind -p N $north_ind -f phased

#####################################################################################

#### I do not know why the above python script provide very high pi values, so I used vcftools to calculate pi, TajD and fst as shown in the following

###Step1.3 vcftools
module load vcftools

vcftools_dir=$vcf_Dir/vcftools

if [ ! -d "$vcftools_dir" ]; then
mkdir -p $vcftools_dir
fi

##Nucleotide diversity
#5kb+2.5kb
#vcftools --gzvcf $Inputvcf --window-pi 5000 --window-pi-step 2500 --keep $vcf_Dir/north.inds --out $vcftools_dir/chr10.N
#vcftools --gzvcf $Inputvcf --window-pi 5000 --window-pi-step 2500 --keep $vcf_Dir/mid.inds --out $vcftools_dir/chr10.M
#vcftools --gzvcf $Inputvcf --window-pi 5000 --window-pi-step 2500 --keep $vcf_Dir/south.inds --out $vcftools_dir/chr10.S
#10kb+5kb
#vcftools --gzvcf $Inputvcf --window-pi 10000 --window-pi-step 5000 --keep $vcf_Dir/north.inds --out $vcftools_dir/chr10.N.10_5kb
#vcftools --gzvcf $Inputvcf --window-pi 10000 --window-pi-step 5000 --keep $vcf_Dir/mid.inds --out $vcftools_dir/chr10.M.10_5kb
#vcftools --gzvcf $Inputvcf --window-pi 10000 --window-pi-step 5000 --keep $vcf_Dir/south.inds --out $vcftools_dir/chr10.S.10_5kb
vcftools --gzvcf $Inputvcf --window-pi 10000 --window-pi-step 5000 --out $vcftools_dir/chr10.all.10_5kb

##Fst
#5kb+2.5kb
#vcftools --gzvcf $Inputvcf --weir-fst-pop $vcf_Dir/north.inds --weir-fst-pop $vcf_Dir/mid.inds --fst-window-size 5000 --fst-window-step 2500 --out $vcftools_dir/chr10_NvsM
#vcftools --gzvcf $Inputvcf --weir-fst-pop $vcf_Dir/south.inds --weir-fst-pop $vcf_Dir/mid.inds --fst-window-size 5000 --fst-window-step 2500 --out $vcftools_dir/chr10_SvsM
#vcftools --gzvcf $Inputvcf --weir-fst-pop $vcf_Dir/north.inds --weir-fst-pop $vcf_Dir/south.inds --fst-window-size 5000 --fst-window-step 2500 --out $vcftools_dir/chr10_NvsS

#10kb+5kb
##vcftools --gzvcf $Inputvcf --weir-fst-pop $vcf_Dir/north.inds --weir-fst-pop $vcf_Dir/mid.inds --fst-window-size 10000 --fst-window-step 5000 --out $vcftools_dir/chr10_NvsM.10_5kb
#vcftools --gzvcf $Inputvcf --weir-fst-pop $vcf_Dir/south.inds --weir-fst-pop $vcf_Dir/mid.inds --fst-window-size 10000 --fst-window-step 5000 --out $vcftools_dir/chr10_SvsM.10_5kb
#vcftools --gzvcf $Inputvcf --weir-fst-pop $vcf_Dir/north.inds --weir-fst-pop $vcf_Dir/south.inds --fst-window-size 10000 --fst-window-step 5000 --out $vcftools_dir/chr10_NvsS.10_5kb



#######################################################################################
###2. Using selscan to estimate the iHS, nSL,XP-EHH
###(An alternative method is to use R-package: REHH2 to estimate the relevant parameters)

module load selscan
module load plink

selscan_dir=$vcf_Dir/selscan

if [ ! -d "$selscan_dir" ]; then
mkdir -p $selscan_dir
fi

tped="/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/plink/SwAsp.chr10.7scaffolds.recode.Potri.rm_het.gt.beagle.tped"
tfam="/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/chr10/vcf_rephase/vcf/rm_het_hwe/beagle/plink/SwAsp.chr10.7scaffolds.recode.Potri.rm_het.gt.beagle.tfam"

ped_dir=`dirname $tped`
#####2.1 Using plink to extract hap data from northern, middle and southern populations separately

#plink --tped $tped --tfam $tfam --keep $ped_dir/north.ind --out $ped_dir/SwAsp.chr10.7scaffolds.north --allow-no-sex --allow-extra-chr --recode transpose
#plink --tped $tped --tfam $tfam --keep $ped_dir/south.ind --out $ped_dir/SwAsp.chr10.7scaffolds.south --allow-no-sex --allow-extra-chr --recode transpose
#plink --tped $tped --tfam $tfam --keep $ped_dir/mid.ind --out $ped_dir/SwAsp.chr10.7scaffolds.mid --allow-no-sex --allow-extra-chr --recode transpose

###2.2 Using perl script to transform the tped file to those files with the information of ancestral state

#sed 's/ /\t/g' $ped_dir/SwAsp.chr10.7scaffolds.north.tped > $ped_dir/temp && mv $ped_dir/temp $ped_dir/SwAsp.chr10.7scaffolds.north.tped
#sed 's/ /\t/g' $ped_dir/SwAsp.chr10.7scaffolds.south.tped > $ped_dir/temp && mv $ped_dir/temp $ped_dir/SwAsp.chr10.7scaffolds.south.tped
#sed 's/ /\t/g' $ped_dir/SwAsp.chr10.7scaffolds.mid.tped > $ped_dir/temp && mv $ped_dir/temp $ped_dir/SwAsp.chr10.7scaffolds.mid.tped

tped_anc="/proj/b2011141/pipeline/perl/SwAsp/tped_anc.01.pl"
anc_Dir="/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/maf005/pos/ancestral"

#perl $tped_anc $ped_dir/SwAsp.chr10.7scaffolds.north.tped  $anc_Dir/SwAsp.pos.outgroup.final.anc
#perl $tped_anc $ped_dir/SwAsp.chr10.7scaffolds.south.tped  $anc_Dir/SwAsp.pos.outgroup.final.anc
#perl $tped_anc $ped_dir/SwAsp.chr10.7scaffolds.mid.tped  $anc_Dir/SwAsp.pos.outgroup.final.anc


###2.3 Using selscan to calculate iHS, nSL, and XP-EHH for northern, southern, and middle populations
###########################################################################################################
###########################################################################################################

##Adding interactive mode to submit the jobs separately and get the job done fast

sel_test=$1
group1=$2  ###north, mid, south or iHS and nSL tests
group2=$3

if [ "$sel_test" == "iHS" ]; then
####iHS
if [ "$group1" == "all" ]; then 
selscan --ihs --tped $ped_dir/SwAsp.chr10.7scaffolds.recode.Potri.rm_het.gt.beagle.ihs.tped --out $selscan_dir/SwAsp.chr10.7scaffolds.$group1 --cutoff 0.05 --threads 2
else
selscan --ihs --tped $ped_dir/SwAsp.chr10.7scaffolds.$group1.ihs.tped --out $selscan_dir/SwAsp.chr10.7scaffolds.$group1 --cutoff 0.05 --threads 2
fi
fi

if [ "$sel_test" == "nSL" ]; then
####nSL
if [ "$group1" == "all" ]; then
selscan --nsl --tped $ped_dir/SwAsp.chr10.7scaffolds.recode.Potri.rm_het.gt.beagle.ihs.tped --out $selscan_dir/SwAsp.chr10.7scaffolds.$group1 --threads 2
else
selscan --nsl --tped $ped_dir/SwAsp.chr10.7scaffolds.$group1.ihs.tped --out $selscan_dir/SwAsp.chr10.7scaffolds.$group1 --threads 2
fi
fi

###XP-EHH  
##Because there need a reference haplotypes file and a test haplotype files, I made three pairs of comparisons N vs. S+M; M vs. S+N; S vs. M+N

###extract .tped files for the two groups of populations
#plink --tped $tped --tfam $tfam --keep $ped_dir/mid_north.ind --out $ped_dir/SwAsp.chr10.7scaffolds.mid_north --allow-no-sex --allow-extra-chr --recode transpose
#plink --tped $tped --tfam $tfam --keep $ped_dir/south_mid.ind --out $ped_dir/SwAsp.chr10.7scaffolds.south_mid --allow-no-sex --allow-extra-chr --recode transpose
#plink --tped $tped --tfam $tfam --keep $ped_dir/south_north.ind --out $ped_dir/SwAsp.chr10.7scaffolds.south_north --allow-no-sex --allow-extra-chr --recode transpose

### Using perl script to transform the tped file to those files with the information of ancestral state

#sed 's/ /\t/g' $ped_dir/SwAsp.chr10.7scaffolds.mid_north.tped > $ped_dir/temp && mv $ped_dir/temp $ped_dir/SwAsp.chr10.7scaffolds.mid_north.tped
#sed 's/ /\t/g' $ped_dir/SwAsp.chr10.7scaffolds.south_mid.tped > $ped_dir/temp && mv $ped_dir/temp $ped_dir/SwAsp.chr10.7scaffolds.south_mid.tped
#sed 's/ /\t/g' $ped_dir/SwAsp.chr10.7scaffolds.south_north.tped > $ped_dir/temp && mv $ped_dir/temp $ped_dir/SwAsp.chr10.7scaffolds.south_north.tped

tped_anc="/proj/b2011141/pipeline/perl/SwAsp/tped_anc.01.pl"
anc_Dir="/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/maf005/pos/ancestral"

#perl $tped_anc $ped_dir/SwAsp.chr10.7scaffolds.mid_north.tped  $anc_Dir/SwAsp.pos.outgroup.final.anc
#perl $tped_anc $ped_dir/SwAsp.chr10.7scaffolds.south_mid.tped  $anc_Dir/SwAsp.pos.outgroup.final.anc
#perl $tped_anc $ped_dir/SwAsp.chr10.7scaffolds.south_north.tped  $anc_Dir/SwAsp.pos.outgroup.final.anc

#if [ "$sel_test" == "XPEHH" ]; then
###the group1 could be: north, mid; the group2 could be south, mid
#north_vs_south
#north_vs_mid
#mid_vs_south
#selscan --xpehh --tped $ped_dir/SwAsp.chr10.7scaffolds.$group1.ihs.tped --tped-ref $ped_dir/SwAsp.chr10.7scaffolds.$group2.ihs.tped --out $selscan_dir/SwAsp.chr10.7scaffolds.${group1}_vs_${group2} --cutoff 0.1  --threads 6
#fi


################################################################################################################################
###3. Calculate H12 for the selected region on chr10
H12_H2H1="/proj/b2011141/tools/SelectionHapStats/scripts/H12_H2H1.py"  ###script used to calculate H12, and H2/H1, it was best used for SNP-based windows instead of physical-based windows, which is because the physical-based windows maybe influenced by local mutation and recombination rate, as discussed in the original paper. Windows with high H12 values are considered as region under hard selective sweeps, whereas windows with both high H12 and high H2/H1 values are considered as region under soft selective sweeps
H12peakFinder="/proj/b2011141/tools/SelectionHapStats/scripts/H12peakFinder.py"
H12_viz="/proj/b2011141/tools/SelectionHapStats/scripts/H12_viz.R"
hapSpectrum="/proj/b2011141/tools/SelectionHapStats/scripts/hapSpectrum_viz.R"
hapData_viz="/proj/b2011141/tools/SelectionHapStats/scripts/hapData_viz.R"
visualize_genome="/proj/b2011141/tools/SelectionHapStats/scripts/visualizeGenomicData.sh"

if [ "$sel_test" == "H12" ]; then

h12_dir=$vcf_Dir/H12
#
if [ ! -d "$h12_dir" ]; then
mkdir -p $h12_dir
fi

if [ "$group1" == "all" ]; then
cut -f 4- $ped_dir/SwAsp.chr10.7scaffolds.recode.Potri.rm_het.gt.beagle.ihs.tped | sed 's/\t/,/g' > $h12_dir/SwAsp.chr10.7scaffolds.$group1.tped  
n_ind=$(cat $ped_dir/SwAsp.chr10.7scaffolds.recode.Potri.rm_het.gt.beagle.tfam |wc -l)
n_hap=$(echo "2*$n_ind" | bc)
python $H12_H2H1 $h12_dir/SwAsp.chr10.7scaffolds.$group1.tped $n_hap -o $h12_dir/SwAsp.chr10.7scaffolds.h12.$group1.200snps.txt -w 200 -j 1 -d 0
else
cut -f 4- $ped_dir/SwAsp.chr10.7scaffolds.$group1.ihs.tped | sed 's/\t/,/g' > $h12_dir/SwAsp.chr10.7scaffolds.$group1.tped
n_ind=$(cat $ped_dir/SwAsp.chr10.7scaffolds.$group1.tfam |wc -l)
n_hap=$(echo "2*$n_ind" |bc)
python $H12_H2H1 $h12_dir/SwAsp.chr10.7scaffolds.$group1.tped $n_hap -o $h12_dir/SwAsp.chr10.7scaffolds.h12.$group1.200snps.txt -w 200 -j 1 -d 0
fi
fi


