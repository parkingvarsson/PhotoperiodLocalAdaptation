#! /bin/bash -l

module load bioinfo-tools java

VCFDir=`dirname $1`
eff_vcf=$1
vcf2maf="/proj/b2011141/pipeline/perl/vcf2maf-master/vcf2maf.pl"
#vcf2maf="/proj/b2011141/tools/vcf2maf-master/vcf2maf.pl"
bgzip="/proj/b2011141/tools/bgzip"

maf_dir=$VCFDir/vcf2maf

if [ ! -d "$maf_dir" ]; then
mkdir -p $maf_dir
fi

out=${eff_vcf##*/}
Outfix=${out%.eff.vcf}

perl $vcf2maf --input-snpeff $eff_vcf --output-maf $maf_dir/$Outfix.eff2.maf && $bgzip $eff_vcf
#perl $vcf2maf --input-vcf $eff_vcf --output-maf $maf_dir/$Outfix.eff2.maf && $bgzip $eff_vcf
#zcat $eff_vcf | perl $vcf2maf --input-snpeff - --output-maf $maf_dir/$Outfix.eff.maf

