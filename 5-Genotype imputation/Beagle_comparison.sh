#! /bin/bash -l



# This job is create several missing genotypes for simulation #

module load perl

# output the file name
#miss_file="/proj/b2011141/nobackup/biyue/snp_asp201/simu_${k}_miss_first100s.vcf"
#packed_file="/proj/b2011141/nobackup/biyue/snp_asp201/Beagle_simu/beagle_simu${k}_imputed_first100s.vcf.gz"
#imputed_file="/proj/b2011141/nobackup/biyue/snp_asp201/Beagle_simu/beagle_simu${k}_imputed_first100s.vcf"

# open several .ped files that has been replace to missing, convert these to .vcf file

for k in {0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5}
do
/proj/b2011141/tools/plink_1.9/\./plink --ped /proj/b2011141/nobackup/biyue/snp_asp201/simu_${k}_miss_first100s.ped --map miss_first100s.map --allow-no-sex --allow-extra-chr --recode vcf --out /proj/b2011141/nobackup/biyue/snp_asp201/simu_${k}_miss_first100s
done

# delete the ./log and .nosex files
rm simu_*_miss_first100s.nosex
rm simu_*_miss_first100s.log

#mkdir Beagle_simu  # make a dir
# do Beagle imputation

for k in {0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5}
do
java -Xmx8000m -jar /proj/b2011141/tools/beagle_4.0/beagle.r1398.jar gtgl=/proj/b2011141/nobackup/biyue/snp_asp201/simu_${k}_miss_first100s.vcf out=/proj/b2011141/nobackup/biyue/snp_asp201/Beagle_simu/beagle_simu${k}_imputed_first100s phase-its=100 impute-its=10
done


# calculate the imputation accuracy by perl
for k in {0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5}
do
zcat /proj/b2011141/nobackup/biyue/snp_asp201/Beagle_simu/beagle_simu${k}_imputed_first100s.vcf.gz >/proj/b2011141/nobackup/biyue/snp_asp201/Beagle_simu/beagle_simu${k}_imputed_first100s.vcf #decompression .gz files and output the same name file
rm /proj/b2011141/nobackup/biyue/snp_asp201/Beagle_simu/beagle_simu${k}_imputed_first100s.vcf.gz # remove .gz files
perl impute_accuracy_beagle_vcf.pl /proj/b2011141/nobackup/biyue/snp_asp201/simu_${k}_miss_first100s.vcf /proj/b2011141/nobackup/biyue/snp_asp201/Beagle_simu/beagle_simu${k}_imputed_first100s.vcf $k #output the calculation data
done


