#! /bin/bash -l

# Mark duplicates in a BAM file


module load bioinfo-tools java
module load samtools/0.1.19

#dedup="/proj/b2010014/GenomePaper/population_genetics/trichocarpa/bwa-mem/realignment/deduplication"
dedup="/proj/b2011141/nobackup/alignments/bwa-mem/asp201/tremula/deduplication"

if [ ! -d "$dedup" ]; then
mkdir -p $dedup
fi

#ln -sf "/home/douglas/tools/picard-tools-1.89/MarkDuplicates.jar" .
MarkDuplicates="/proj/b2011141/tools/picard-tools-1.115/MarkDuplicates.jar"

if [ "$1" == "" -o "$2" != "" ] ; then
	echo One BAM file at a time
	exit 1
fi

InputBAM=$1

BAM_index=${1%.bam}.bai
if [ -f $BAM_index ] ; then
  rm -f $BAM_index
fi
(samtools index $1 ) &
echo waiting for indexing of BAM files
wait
echo indexing completed

BAMs=${1##*/}
OutputBAM=$dedup/${BAMs%.bam}.dedup.bam
Metrics=$dedup/${BAMs%.bam}.metrics

# java -Xmx?
java -Xmx2500m -jar $MarkDuplicates INPUT=$InputBAM OUTPUT=$OutputBAM METRICS_FILE=$Metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

echo deduplication completed                                    


