##Simulate selective sweeps using msms and analöyse using SweepFinder2 to assess region affected by sweep
#! /bin/bash -l


/proj/b2011141/tools/msms/bin/msms -ms 52 1000 -t 5670 -r 1298 -SAA 3100 -SaA 1550 -N 92000 -SF 0 -Sp 0.5 | ./ms2sf2.pl

for i in sf2_sim*.sf2
do
    fout=`basename $i ".sf2"`
    /home/pari/src/SF2/SweepFinder2 -sg 2200 $i $fout".out"
done
