#!/bin/bash -l


set -e
set -x

module load bioinfo-tools
module load FastQC/0.10.1

PROJ="/proj/b2011141"
suffixIn=fastq
needsTranslate=yes  # or nothing

TmpID=${SLURM_JOB_ID:-$$}

SCRIPTS="$PROJ/pipeline"
TOOLS="$PROJ/tools"
cutadapt="$TOOLS/cutadapt"
trimmomatic="$TOOLS/trimmomatic-0.30.jar"

CutadaptConf="$SCRIPTS/Cutadapt_final.conf"
AdapterFasta="$SCRIPTS/adapter.fa"

fileIn=${1##*/}                     # SwAsp006_650bp_1.fastq.gz
filePath=${1%/*}                    # /proj/b2011141/sequenceData/tremula/SwAsp006
genericIn=${1%_1.${suffixIn}.gz}    # /proj/b2011141/sequenceData/tremula/SwAsp006/SwAsp006_650bp
speciesP=${1##*sequenceData/}       # tremula/SwAsp006/SwAsp006_650bp_1.fastq.gz
species=${speciesP%%/*}             # tremula
sampleFile=${fileIn%_1.${suffixIn}.gz}    # SwAsp006_650bp

echo suffixIn=$suffixIn
echo fileIn=$fileIn
echo filePath=$filePath
echo genericIn=$genericIn
echo speciesP=$speciesP
echo species=$species
echo sampleFile=$sampleFile

if [ ! -e "${genericIn}_1.${suffixIn}.gz" ] ; then
        if [ -L "${genericIn}_1.${suffixIn}.gz" ]; then
                 echo link to read 1 file ${genericIn}_1.${suffixIn}.gz is broken ; exit 1
        else
                 echo could not find read 1 file ${genericIn}_1.${suffixIn}.gz ; exit 1
        fi
    fi
if [ ! -e "${genericIn}_2.${suffixIn}.gz" ] ; then
        if [ -L "${genericIn}_2.${suffixIn}.gz" ]; then
                 echo link to read 2 file ${genericIn}_2.${suffixIn}.gz is broken ; exit 1
        else
                 echo could not find read 1 file ${genericIn}_2.${suffixIn}.gz ; exit 1
        fi
fi

TRIM_DIR="$PROJ/nobackup/sequenceData/QCfiltered/$species"
#TRIM_DIR="$PROJ/nobackup/private/new_nobackup_area/QCfiltered/$species"

if [ ! -d "$TRIM_DIR" ]; then
         mkdir $TRIM_DIR
fi

RawR1=${genericIn}_1.${suffixIn}.gz
RawR2=${genericIn}_2.${suffixIn}.gz

echo read 1 $RawR1
echo read 2 $RawR2



# detect Phred quality scoring

QualityScoringOption=
Quality=$($TOOLS/phredDetector.pl $RawR1)
if [ "$Quality" = "33" ] ; then
        QualityScoringOption="-phred33"
elif [ "$Quality" = "64" ] ; then
       QualityScoringOption="-phred64"
else
echo "phredDetector couldn't autodetect quality for $RawR1, return value was '$Quality'"
exit 1
fi
echo "phredDetector infers $RawR1 needs quality scorring option $QualityScoringOption"


# create a FastQC report for the data

RawFastQCResults="/proj/b2011141/nobackup/sequenceData/fastQC/fastqc-results"
mkdir -p $RawFastQCResults

echo running FastQC in parallel on both files of raw reads, placing results in $RawFastQCResults

fastqc --quiet --threads 2 --outdir $RawFastQCResults --extract $RawR1 $RawR2

# directories containing the extracted files are ${RawR1%.${suffixIn}}_fastqc and ${RawR2%.${suffixIn}}_fastqc
# within each directory, the file Images/per_base_quality.png holds the read quality plots
RawR1_Qual=$RawFastQCResults/${sampleFile}_1_fastqc/Images/per_base_quality.png
RawR2_Qual=$RawFastQCResults/${sampleFile}_2_fastqc/Images/per_base_quality.png

# run trimmomatic

Trimmomatic_OUT_I_1=$TRIM_DIR/$sampleFile.all-trimmomatic.$TmpID.1.fq.gz
Trimmomatic_OUT_SE_1=$TRIM_DIR/$sampleFile.all-trimmomatic.$TmpID.forward.se.fq.gz
Trimmomatic_OUT_I_2=$TRIM_DIR/$sampleFile.all-trimmomatic.$TmpID.2.fq.gz
Trimmomatic_OUT_SE_2=$TRIM_DIR/$sampleFile.all-trimmomatic.$TmpID.reverse.se.fq.gz
Trimmomatic_OUT_SE=$TRIM_DIR/$sampleFile.all-trimmomatic.$TmpID.se.fq.gz

Trimmomatic_OUT_I_1_final=$TRIM_DIR/$sampleFile.all-trimmomatic.1.fq.gz
Trimmomatic_OUT_I_2_final=$TRIM_DIR/$sampleFile.all-trimmomatic.2.fq.gz
Trimmomatic_OUT_SE_final=$TRIM_DIR/$sampleFile.all-trimmomatic.se.fq.gz

# do not remove any of the intermediate files for now

echo "trimmomatic adapter + quality trim, $translate_OUT_1 $translate_OUT_2 ..."
trimmomatic_Output=$TRIM_DIR/$sampleFile.trimOutput
echo "trimmomatic output to $trimmomatic_Output"

java -Xmx20g -jar $trimmomatic PE -threads 8 $QualityScoringOption $RawR1 $RawR2 $Trimmomatic_OUT_I_1 $Trimmomatic_OUT_SE_1 $Trimmomatic_OUT_I_2 $Trimmomatic_OUT_SE_2 ILLUMINACLIP:${AdapterFasta}:1:30:9 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36 2>&1 | tee $trimmomatic_Output 1>&2

zcat $Trimmomatic_OUT_SE_1 $Trimmomatic_OUT_SE_2 | gzip -c - > $Trimmomatic_OUT_SE && rm -f $Trimmomatic_OUT_SE_1 $Trimmomatic_OUT_SE_2
rm -f $translate_OUT_1 $translate_OUT_2

mv -f $Trimmomatic_OUT_I_1 $Trimmomatic_OUT_I_1_final
mv -f $Trimmomatic_OUT_I_2 $Trimmomatic_OUT_I_2_final
mv -f $Trimmomatic_OUT_SE  $Trimmomatic_OUT_SE_final

# learn what Trimmomatic did
Stats=$(grep "^Input Read Pairs:" $trimmomatic_Output)
NInput=$(echo $Stats | cut -d" " -f4)
NPairSurviving=$(echo $Stats | cut -d" " -f7-8)
NForwardSurviving=$(echo $Stats | cut -d" " -f12-13)
NReverseSurviving=$(echo $Stats | cut -d" " -f17-18)
NDropped=$(echo $Stats | cut -d" " -f20-21)


# Now do a FastQC report on the post-QC reads

TrimmedFastQCResults="$TRIM_DIR/fastqc-results"
mkdir -p $TrimmedFastQCResults

echo running FastQC in parallel on both files of post-QC reads, placing results in $TrimmedFastQCResults

sampleTrimFile=$sampleFile.all-trimmomatic

fastqc --quiet --threads 3 --outdir $TrimmedFastQCResults --extract $Trimmomatic_OUT_I_1_final $Trimmomatic_OUT_I_2_final $Trimmomatic_OUT_SE_final

# directories containing the extracted files are ${RawR1%.gz}_fastqc and ${RawR2%.gz}_fastqc
# within each directory, the file Images/per_base_quality.png holds the read quality plots
TrimmedR1_Qual=$TrimmedFastQCResults/$sampleTrimFile.1.fq_fastqc/Images/per_base_quality.png
TrimmedR2_Qual=$TrimmedFastQCResults/$sampleTrimFile.2.fq_fastqc/Images/per_base_quality.png
TrimmedSE_Qual=$TrimmedFastQCResults/$sampleTrimFile.se.fq_fastqc/Images/per_base_quality.png

echo Producing HTML file with FastQC per-base-quality results

HTMLResults="$TrimmedFastQCResults/HTML"
mkdir -p $HTMLResults

RawR1_Qual_final=$HTMLResults/$sampleFile.1.per_base_quality.png
RawR2_Qual_final=$HTMLResults/$sampleFile.2.per_base_quality.png
TrimmedR1_Qual_final=$HTMLResults/$sampleTrimFile.1.per_base_quality.png
TrimmedR2_Qual_final=$HTMLResults/$sampleTrimFile.2.per_base_quality.png
TrimmedSE_Qual_final=$HTMLResults/$sampleTrimFile.se.per_base_quality.png
Qual_final=$sampleFile.complete.per_base_quality.png

convert -resize 350 $RawR1_Qual $RawR1_Qual_final
convert -resize 350 $RawR2_Qual $RawR2_Qual_final
convert -resize 350 $TrimmedR1_Qual $TrimmedR1_Qual_final
convert -resize 350 $TrimmedR2_Qual $TrimmedR2_Qual_final
convert -resize 350 $TrimmedSE_Qual $TrimmedSE_Qual_final
convert \( $RawR1_Qual_final $RawR2_Qual_final \
           \( -background white -pointsize 14 label:"$sampleFile\n\nQuality scoring option: $QualityScoringOption\n\nInput pairs: $NInput\nPairs surviving: $NPairSurviving\nForward surviving: $NForwardSurviving\nReverse surviving: $NReverseSurviving\nDropped: $NDropped\n\n---------\n\nRaw R1           Raw R2\n\nTrimmed R1    Trimmed R2    Trimmed SE" -gravity west \) \
           +append \) \
        \( $TrimmedR1_Qual_final $TrimmedR2_Qual_final $TrimmedSE_Qual_final +append \) \
        -background white -append $HTMLResults/$Qual_final
    rm -f $RawR1_Qual_final $RawR2_Qual_final $TrimmedR1_Qual_final $TrimmedR2_Qual_final $TrimmedSE_Qual_final

    cat <<__report__ > $HTMLResults/$sampleFile.html
    <h3>${sampleFile} FastQC Results</h3>
    <img src="$Qual_final">
    <table border="0" cellspacing="20">
    <tr> <th>Sample</th> <th>Input pairs</th> <th>Pairs surviving</th> <th>Forward surviving</th> <th>Reverse surviving</th> <th>Dropped</th> </tr>
    <tr> <td>$sampleFile</td> <td>$NInput</td> <td>$NPairSurviving</td> <td>$NForwardSurviving</td> <td>$NReverseSurviving</td> <td>$NDropped</td> </tr>
    </table>
    __report__
