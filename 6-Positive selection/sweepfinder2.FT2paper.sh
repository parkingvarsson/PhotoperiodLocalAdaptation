#! /bin/bash -l


###The main aim of this script is use SweepFinder2 to calculate CLR on the ~700 region on Chr10 and across the genome, since this script is based on Derived alllele frequency across the whole genome, so first I need to do is to extract all polymorhpic sites from the seven scaffolds in the FT2 region

module load bioinfo-tools
module load plink

###the inputped file is the tped file with information of ancestral (0) and derived (1) and also contain all polymorhpic sites (7734522 in total)
Inputped="/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/all/beagle/plink/SwAsp_94samples.filter.gt.beagle.ihs.tped"
Inputfam="/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/vcf/all/beagle/plink/SwAsp_94samples.filter.gt.beagle.ihs.tfam"

Out=${Inputped##*/}
echo $Out
Out_index=${Out%%.tped}

ped_Dir=`dirname $Inputped`

#######################################################################################
###1. use plink to extract and create tped files seprately for north, mid and south populations

module load plink

group_dir=$ped_Dir/group

if [ ! -d "$group_dir" ]; then
mkdir -p $group_dir
fi

#plink --tped $Inputped --tfam $Inputfam --keep $ped_Dir/north.ind --out $group_dir/${Out_index}.north --allow-no-sex --allow-extra-chr --recode transpose --missing-genotype 'O'
#plink --tped $Inputped --tfam $Inputfam --keep $ped_Dir/mid.ind --out $group_dir/${Out_index}.mid --allow-no-sex --allow-extra-chr --recode transpose --missing-genotype 'O'
#plink --tped $Inputped --tfam $Inputfam --keep $ped_Dir/south.ind --out $group_dir/${Out_index}.south --allow-no-sex --allow-extra-chr --recode transpose --missing-genotype 'O'


#######################################################################################
###2. Use R script to sumarize the derived allele account for each SNP and also extract the 7 scaffolds from chr10 region as an independent dataset and also exclude those 7 scaffolds from the genome-wide datasets

sf2_r="/proj/b2011141/pipeline/R/local_adaptation_paper/sweepfinder2/chr10.sweepfinder2.R"   ###the script creating input of SweepFinder2 from tped file
sf2_chr_r="/proj/b2011141/pipeline/R/local_adaptation_paper/sweepfinder2/genome.chr.sweepfinder2.R"   ###the script creating input of SweepFinder2 from tped file for each chromosome
sf_dir="/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/SweepFinder2"
sf_chr_dir="/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/SweepFinder2/chromosome"

if [ ! -d "$sf_chr_dir" ]; then
mkdir -p $sf_chr_dir
fi

#Rscript $sf2_r $group_dir ${Out_index}.north.tped $sf_dir 
#Rscript $sf2_r $group_dir ${Out_index}.mid.tped $sf_dir 
#Rscript $sf2_r $group_dir ${Out_index}.south.tped $sf_dir 

#Rscript $sf2_chr_r $group_dir ${Out_index}.north.tped $sf_chr_dir
#Rscript $sf2_chr_r $group_dir ${Out_index}.mid.tped $sf_chr_dir
#Rscript $sf2_chr_r $group_dir ${Out_index}.south.tped $sf_chr_dir



#######################################################################################
###3. Using SweepFinder to run 

sf2="/proj/b2011141/tools/SF2/SweepFinder2"

###3.1 Compute empirical frequency spectrum (-f) for different population groups separately
#$sf2 -f $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.north.tped.genome.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.north.spect.sf2
#$sf2 -f $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.mid.tped.genome.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.mid.spect.sf2
#$sf2 -f $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.south.tped.genome.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.south.spect.sf2

###3.2 scan for selective sweeps
###10kb window
#$sf2 -lg 10000 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.north.tped.chr10.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.north.spect.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.north.chr10.10kb.sf2.out
#$sf2 -lg 10000 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.mid.tped.chr10.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.mid.spect.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.mid.chr10.10kb.sf2.out
#$sf2 -lg 10000 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.south.tped.chr10.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.south.spect.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.south.chr10.10kb.sf2.out

#2kb
#$sf2 -lg 2000 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.north.tped.chr10.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.north.spect.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.north.chr10.2kb.sf2.out
#$sf2 -lg 2000 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.mid.tped.chr10.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.mid.spect.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.mid.chr10.2kb.sf2.out
#$sf2 -lg 2000 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.south.tped.chr10.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.south.spect.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.south.chr10.2kb.sf2.out

###3.3 scan for selective sweep for each chromosome

#2kb
chr=$1 ##01..19
pop=$2 ##mid, north, south
$sf2 -lg 2000 $sf_chr_dir/SwAsp_94samples.filter.gt.beagle.ihs.$pop.tped.PseudoChr$chr.chr_tped.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.$pop.spect.sf2 $sf_chr_dir/SwAsp_94samples.filter.gt.beagle.ihs.$pop.tped.PseudoChr$chr.chr_tped.sf2.out
#$sf2 -lg 2000 $sf_chr_dir/SwAsp_94samples.filter.gt.beagle.ihs.north.tped.PseudoChr$chr.chr_tped.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.north.spect.sf2 $sf_chr_dir/SwAsp_94samples.filter.gt.beagle.ihs.north.tped.PseudoChr$chr.chr_tped.sf2.out
#$sf2 -lg 2000 $sf_chr_dir/SwAsp_94samples.filter.gt.beagle.ihs.south.tped.PseudoChr$chr.chr_tped.sf2 $sf_dir/SwAsp_94samples.filter.gt.beagle.ihs.north.spect.sf2 $sf_chr_dir/SwAsp_94samples.filter.gt.beagle.ihs.south.tped.PseudoChr$chr.chr_tped.sf2.out




