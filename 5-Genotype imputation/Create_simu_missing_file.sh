#! /bin/bash -l

#set -e

# This job is for simulating by Beagle and calculate accuracy#

#module load bioinfo-tools
#module load samtools/0.1.19
#module load bwa

ori_file="/proj/b2011141/nobackup/biyue/snp_asp201/miss_first100s.ped"
out_file="/proj/b2011141/nobackup/biyue/snp_asp201/miss_first100s_re.ped"
perl_script="/proj/b2011141/pipeline/biyue/renew_pedForm.pl"

perl $perl_script $ori_file $out_file #change genotype from 'G G' to 'G/G'

Rscript /proj/b2011141/pipeline/biyue/createFile_imputate.R #create many files with different missing percentages

# replace some characters in the files
for file in simu_*_miss_first100s_re.ped
do
cat $file | sed -e 's/\// /g; s/\"//g' >$file.ped
done

# delect the files generated from R script
rm simu_*_miss_first100s_re.ped

# rename the files to a make sense one
rename _re.ped.ped .ped *.ped.ped


