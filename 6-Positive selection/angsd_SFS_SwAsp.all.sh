#! /bin/bash -l


angsd="/proj/b2011141/tools/angsd/angsd"
emOptim2="/proj/b2011141/tools/angsd0.602/misc/emOptim2"
thetaStat="/proj/b2011141/tools/angsd0.602/misc/thetaStat"
bam_list_north="/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/ANGSD/bam/northern/northern.bam.list"
bam_list_middle="/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/ANGSD/bam/middle/middle.bam.list"
bam_list_south="/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/ANGSD/bam/southern/southern.bam.list"

ref="/proj/b2011141/nobackup/reference/asp201/Potra01-genome.fa"
OutDir="/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/ANGSD/SFS/all"
region="/proj/b2011141/nobackup/all_populations/GATK/HC/bed_file/final_version/Potra/split"
chr="$region/Potra.SwAsp.filter.$1.region"

nInd=$(cat $bam_list_all | wc -l)
#nChrom=$(echo "2*$nInd" | bc)
nChrom=$nInd

echo $nInd
echo $nChrom

if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

#first generate .saf file
$angsd -bam $bam_list_all -minMapQ 30 -minQ 20 -GL 1 -doSaf 1 -out $OutDir/SwAsp_94.$1 -anc $ref -fold 1 -rf $chr 

#use emOptim2 to optimization
#$emOptim2 $OutDir/SwAsp_94.$1.saf $nChrom -maxIter 100 -P 4 > $OutDir/SwAsp_94.$1.sfs

#calculate thetas
#$angsd -bam $bam_list_all -out $OutDir/SwAsp_94.$1 -doThetas 1 -GL 1 -doSaf 1 -anc $ref -rf $chr -pest $OutDir/SwAsp_94.$1.sfs -minMapQ 30 -minQ 20 -fold 1 

#calculate Tajimas
#$thetaStat make_bed $OutDir/SwAsp_94.$1.thetas.gz
#$thetaStat do_stat $OutDir/SwAsp_94.$1.thetas.gz -nChr $nChrom
$thetaStat do_stat $OutDir/SwAsp_94.$1.thetas.gz -nChr $nChrom -win 1000 -step 1000 -outnames $OutDir/SwAsp_94.$1.thetas1kbwindow.gz
#$thetaStat do_stat $OutDir/SwAsp_94.$1.thetas.gz -nChr $nChrom -win 5000 -step 5000 -outnames $OutDir/SwAsp_94.$1.thetas5kbwindow.gz
#$thetaStat do_stat $OutDir/tremula_$1.thetas.gz -nChr $nChrom -win 10000 -step 10000 -outnames $OutDir/tremula_$1.thetas10kbwindow.gz
#$thetaStat do_stat $OutDir/tremula_$1.thetas.gz -nChr $nChrom -win 100000 -step 100000 -outnames $OutDir/tremula_$1.thetas100kbwindow.gz
#$thetaStat do_stat $OutDir/tremula_$1.thetas.gz -nChr $nChrom -win 500000 -step 500000 -outnames $OutDir/tremula_$1.thetas500kbwindow.gz
#$thetaStat do_stat $OutDir/tremula_$1.thetas.gz -nChr $nChrom -win 10000 -step 5000 -outnames $OutDir/tremula_$1.thetas10kbwindow.5kbsliding.gz
#$thetaStat do_stat $OutDir/tremula_$1.thetas.gz -nChr $nChrom -win 100000 -step 10000 -outnames $OutDir/tremula_$1.thetas100kbwindow.10kbsliding.gz
#$thetaStat do_stat $OutDir/tremula_$1.thetas.gz -nChr $nChrom -win 500000 -step 50000 -outnames $OutDir/tremula_$1.thetas500kbwindow.50kbsliding.gz

#calculate MAF

#calculate posterior probabilities of sample allele frequencies
#$angsd -bam $bam_list_tremula -GL 1 -doSaf 1 -anc $ref -rf $chr -minMapQ 30 -minQ 20 -pest $OutDir/tremula_$1.sfs -out $OutDir/tremula_$1.rf


