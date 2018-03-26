#! /usr/bin/perl -w
### The aim of this script is to calculate the imputation accuracy by each sample and each SNP ###
# input patten are (1) input missing file name; (2) input imputed file name; (3) output pattern, like "0.05"
use strict;
use Statistics::Basic qw(:all);

#open ORIGINAL, "<", $ARGV[0] or die $!;   # open a file end with .ped, usually, this file is "BD_chr6.ped"
#or above line also as:
#open ORIGINAL, "</proj/b2011141/nobackup/biyue/snp_asp201/miss_first100s.vcf" or die $!;
open ORIGINAL, "miss_first100s.vcf" or die $!;
open MISS, "<", $ARGV[0] or die $!;  # open another file that has been random missing setted, In this case, file as "simu_*_miss_first100s.vcf"
open TEST, "<", $ARGV[1] or die $!;  # open a corresponding imputed file. In this case, file as "beagle_simu*_imputed_first100s.vcf"
#open SAMPLE, ">", join("", "/proj/b2011141/nobackup/biyue/snp_asp201/Beagle_simu/sample_",$ARGV[2], "_error_rate.txt") or die $!;  # create a file which including error rate of each sample
open SAMPLE, ">", join("", $ARGV[2], "_sample_error_rate.txt") or die $!;
#open SNP, ">", join("", "/proj/b2011141/nobackup/biyue/snp_asp201/Beagle_simu/SNP_",$ARGV[2], "_error_rate.txt") or die $!;  # create a file which including error rate of each SNP
open SNP, ">", join("", $ARGV[2], "_SNP_error_rate.txt") or die $!;
open R2, ">", join("", $ARGV[2], "_R2_rate.txt") or die $!;

my @snp_title = ();
my @sample_title=();
my @allele_r2=();
my @dosage_r2=();
my %miss;
my %test;
my %impute;
my %snp_all;
my %snp_false;
my %sample_all;
my %sample_false;
my $row1 = 0;
my $row2 = 0;
my $row3 = 0;
my $miss_all=0;
my $false_all=0;
my $col=0;



while(<MISS>){
    chomp;
    next if (/^##/);
    if (/^#CHROM/){
    	@sample_title=split("\t",$_);
    	@sample_title=@sample_title[9..$#sample_title];
        next;
    }
    $row2++;
    my @line=split("\t", $_);
    my $title=join("_",$line[0],$line[1]);
    push (@snp_title, $title);
    @line=@line[9..$#line];
    my $col1 = 0;
    foreach(@line){ # search which row and column is missing and create a hash:position is key and ./. is value
        $col1++;
        next if ($_ ne "./.");
        my $name = join("_", $row2, $col1);
        $miss{$name}=$_;  
    }
}

while(<ORIGINAL>){
		chomp;
		next if (/^#/);
		$row1++;
		my @line=split("\t", $_);
		my $ref=$line[3];
		my $alt=$line[4];
		@line=@line[9..$#line];
		my $col1 = 0;
		foreach(@line){ # search which SNP is artificially delete and recode the true value
            $col1++;
            next if ($_ eq "./.");
            my $name = join("_", $row1, $col1);
            if (exists $miss{$name}){
            	my @allele=split("/",$_);
            	foreach(@allele) {
            		if ($_ eq "0") {
            			$_=$ref;
            		} elsif ($_ eq "1"){
            			$_=$alt;
            		}
            	}
            	@allele=sort @allele;
            	my $new_allele=join("/",@allele);
                $impute{$name}=$new_allele;
            }
        }
    $col=$col1;
}
my @imputed=keys %impute;


#print IMPUTE join(" ", "Row","Colume","SNP"), "\n";
#foreach (keys %impute){
#		my @name_p = split("_", $_);
#		print IMPUTE join(" ", @name_p, $impute{$_}), "\n";
#}
print R2 join("\t", "SNP","Allelic R2","Dosage R2");

while(<TEST>){
	chomp;
	next if (/^#/);
	$row3++;
    $_=~ s/\|/\//g;
    #print $_, "\n";
    my @line=split("\t", $_);
    my $title=join("_",$line[0],$line[1]);
    my $ref=$line[3];
    my $alt=$line[4];
    my @r2=split(";",$line[7]);
    $r2[0]=~ s/AR2=//;
    $r2[1]=~ s/DR2=//;
    push(@allele_r2, $r2[0]);
    push(@dosage_r2, $r2[1]);
    print R2 join("\t", $title,$r2[0],$r2[1]), "\n";
    @line=@line[9..$#line];
    my $col1 = 0;
    foreach (@line){
        $col1++;
        my $name = join("_", $row3, $col1);
        if (exists $impute{$name}) { # search inputed artificially deleted SNPs, and recode the imputed value
            my @info=split(":", $_);
            my @allele=split("/", $info[0]);
            foreach(@allele){
                if ($_ eq "0") {
                    $_=$ref;
                } elsif ($_ eq "1"){
                    $_=$alt;
                }
            }
            @allele=sort @allele;
            my $new_allele=join("/",@allele);
            $test{$name}=$new_allele;
        }
    }
}

my $mean_allele=mean(@allele_r2);
my $mean_dosage=mean(@dosage_r2);

print R2 join("\t","Mean", $mean_allele, $mean_dosage), "\n";

#print IMPUTE join(" ", "Row","Colume","SNP"), "\n";
#foreach (keys %test){
#		my @name_p = split("_", $_);
#		print IMPUTE join(" ", @name_p, $test{$_}), "\n";
#}

#print FALSE join("\t", "SNP","Sample","SNP_true","SNP_imputed"), "\n";
foreach (@imputed){
	next if ($test{$_} eq $impute{$_});
	my @name_p = split("_", $_);
	#print FALSE join("\t", @name_p, $impute{$_}, $test{$_}), "\n";
   $false_all++;
}
#print $false_all, "\n";

$miss_all=@imputed;
#print $miss_all, "\n";
my $rate_all = $false_all/$miss_all;

foreach(@imputed){
    my @name_p = split("_", $_);
    my $snp_p=$name_p[0];
    $snp_all{$snp_p}++;
    next if ($test{$_} eq $impute{$_});
    $snp_false{$snp_p}++;
}

foreach(@imputed){
    my @name_p = split("_", $_);
    my $sample_p=$name_p[1];
    $sample_all{$sample_p}++;
    next if ($test{$_} eq $impute{$_});
    $sample_false{$sample_p}++;
}

print SNP join("\t", "SNP","num/false", "num/miss", "ratio"), "\n";
my $num1=0;
while($num1<=$row1){
    $num1++;
    next unless (exists $snp_all{$num1});
    my $rate=0;
    my $snp_false=0;
    if (exists $snp_false{$num1}){
        $rate=$snp_false{$num1}/$snp_all{$num1};
        $snp_false=$snp_false{$num1};
    } else {
        $rate=0/$snp_all{$num1};
    }
    print SNP join("\t", $snp_title[$num1-1],$snp_false,$snp_all{$num1},$rate), "\n";
}
print SNP join("\t", "All", $false_all,$miss_all,$rate_all),"\n";

print SAMPLE join("\t", "Sample","num/false", "num/miss", "ratio"), "\n";
my $num2=0;
while($num2<=$col){
    $num2++;
    next unless (exists $sample_all{$num2});
    my $rate=0;
    my $sample_false=0;
    if (exists $sample_false{$num2}){
        $rate=$sample_false{$num2}/$sample_all{$num2};
        $sample_false=$sample_false{$num2};
    } else {
        $rate=0/$sample_all{$num2};
    }
    print SAMPLE join("\t", $sample_title[$num2-1],$sample_false,$sample_all{$num2},$rate), "\n";
}
print SAMPLE join("\t", "All", $false_all,$miss_all,$rate_all),"\n";


