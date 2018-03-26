#! /bin/bash -l



#module load bioinfo-tools java
#module load samtools
#module load GATK

samtools="/proj/b2011141/tools/samtools-0.1.19/samtools"
REF="/proj/b2011141/nobackup/reference/asp201/Potra01-genome.fa"
AlignmentDir="/proj/b2010014/nobackup/population_genetics/tremula/aspen_snp_calling/bwa-mem/tremula/realignment/deduplication"
OutDir="/proj/b2010014/nobackup/population_genetics/tremula/aspen_snp_calling/GATK/HC/tremula"
bed="/proj/b2011141/nobackup/reference/asp201/bed/split/Potra.$2.bed"
GATK="/proj/b2011141/tools/GATK/3.2.2/GenomeAnalysisTK.jar"

InputBAMs=$AlignmentDir/$1.merge.dedup.bam

InputBAMs_index=${InputBAMs}.bai
if [ ! -f $InputBAMs_index ]; then
      	$samtools index $InputBAMs
fi

echo InputBAMs=\"$InputBAMs\"


Outfile="$1.$2.gatk.hap.vcf"

java -jar $GATK -R $REF \
 -T HaplotypeCaller \
 -I $InputBAMs \
 --emitRefConfidence GVCF \
 --variant_index_type LINEAR \
 --variant_index_parameter 128000 \
 --allow_potentially_misencoded_quality_scores \
 -L $bed \
 --heterozygosity 0.015 \
 --indel_heterozygosity 0.0025 \
 -o $OutDir/$Outfile

#$Tools/bgzip $OutDir/$Outfile
#$Tools/tabix -p vcf $OutDir/$Outfile.gz

echo Haplotype calling completed


