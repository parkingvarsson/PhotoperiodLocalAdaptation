#! /bin/bash -l

module load bioinfo-tools java



VCFDir=`dirname $1`
eff_vcf=`basename $1 .vcf.gz`.eff.vcf
snpEff="/proj/b2011141/tools/snpEff_3.6/snpEff/snpEff.jar"
bgzip="/proj/b2011141/tools/bgzip"
tabix="/proj/b2011141/tools/tabix"

# -fi intervals.bed # only annotate in regions of interval file
# -ud size_in_bases #this option can be used to change the default upstream and downstream interval size
# -t mutlithreaded option # this allows to use several cores available in the local machine but can make some features not available, e.g. statistics
# -v verbose option #this shows all chromosome names and their respective lengths
# -spliceSiteSize size_in_bases


eff=$VCFDir/eff_v1.1

if [ ! -d "$eff" ]; then
mkdir -p $eff
fi


#zcat $1 | java -Xmx8g -jar $snpEff eff -sequenceOntology -v -ud 2000 -spliceSiteSize 4 GAM_asp201-001 - > $eff/$eff_vcf
zcat $1 | java -Xmx8g -jar $snpEff eff -sequenceOntology -v -ud 2000 -spliceSiteSize 4 Potra-v1.1 - > $eff/$eff_vcf
$bgzip $eff/$eff_vcf
$tabix -p vcf $eff/$eff_vcf.gz

