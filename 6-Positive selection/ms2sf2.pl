#!/usr/bin/perl

$S0=0;
$line=<>;
$num=(split/\s+/,$line)[1];
##$sites=(split/\s+/,$line)[3];
$sites=500000;
$line=<>;
$line=<>;
my @segsites;
$filenum=0;

while(<>){
    chomp;
    if (/^\/\//){
	$line=<>;
	chomp $line;
	$S=(split/\s+/,$line)[1];
	$line=<>;
	@pos = split/\s/,$line;

    }elsif (/^$/) {
	$filenum++;
	$fout = "sf2_sim".$filenum.".sf2";
	open(my $fh,'>',$fout) or die("Couuld not open $fout for writing\n");
	print $fh "position\tx\tn\tfolded\n";
	for($i=0;$i<$S;$i++){
	    print $fh int($pos[$i+1]*$sites+0.5),"\t";
	    print $fh $segsites[$i],"\t";
	    print $fh $num,"\t0\n";
	}
	close $fh;
	@segsites=();
    }else{
	@line = split//,$_;
	for($i=0;$i<$S;$i++){
	   $segsites[$i]+=$line[$i];
	}
    }
}
$filenum++;
$fout = "sf2_sim".$filenum.".sf2";
open(my $fh,'>',$fout) or die("Couuld not open $fout for writing\n");
print $fh "position\tx\tn\tfolded\n";
for($i=0;$i<$S;$i++){
    print $fh int($pos[$i+1]*$sites+0.5),"\t";
    print $fh $segsites[$i],"\t";
    print $fh $num,"\t0\n";
}
close $fh;
@segsites=();
