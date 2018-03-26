#! /bin/bash -l

#SBATCH -A b2011141
#SBATCH -p core
#SBATCH -o GEMMA.lmm.noPC.out
#SBATCH -e GEMMA.lmm.noPC.err
#SBATCH -J GEMMA.lmm.noPC.job
#SBATCH -t 6:00:00
#SBATCH --mail-user jing.wang@emg.umu.se
#SBATCH --mail-type=ALL


gemma="/proj/b2011141/tools/gemma"
genotype_folder="/proj/b2011141/nobackup/local_adaptation_paper/asp201_94/budset/GEMMA"

########Step2.kinship matrix estimate by GEMMA

### Calculate kinship matrix using GEMMA, and using centered relatedness matrix function
cd $genotype_folder
$gemma -bfile SwAsp.94samples.beagle.maf.recode_part1 -gk 1 -o gemma_centered_kinship
 
########Step 3.do GWAS with GEMMA #####

##Second, start to use GEMMA to do GWAS for each phenotype
##because we divided the input bed file into 4 parts, so I use the for loop to do GWAS sequentially for each file
kinship=$genotype_folder/output/gemma_centered_kinship.cXX.txt

cd $genotype_folder
for k in {1..4}
do
##lmm 4 here means we perfer to use three methods to do frequentist analysis choice: 1.wald test, 2. likelihood ratio test, 3. score test
$gemma -bfile $genotype_folder/SwAsp.94samples.beagle.maf.recode_part${k} -k $kinship -n 1 -lmm 4 -o SwAsp.94.filter_part${k}.budset.noPC 
done


###combine the 4 part of each traits into one file
for file2 in {2..4}
do
sed '1d' $genotype_folder/output/SwAsp.94.filter_part$file2.budset.noPC.assoc.txt > temp && mv temp $genotype_folder/output/SwAsp.94.filter_part$file2.budset.noPC.assoc.txt
done

input=
for file2 in {1..4}
do
input="$input $genotype_folder/output/SwAsp.94.filter_part$file2.budset.noPC.assoc.txt"
done
cat $input > $genotype_folder/output/SwAsp.94.filter.all.budset.noPC.assoc.txt 
rm $input


