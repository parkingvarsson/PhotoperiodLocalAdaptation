Title: Ekebo_Savar_budset.txt.

#site	year	pos	block	row	plant	clone	budset
#Savar	2006	1	1	1	1	46	205
#Savar	2006	2	1	1	2	96	NA
#Savar	2006	3	1	1	3	28	226
#Savar	2006	4	1	1	4	83	226

 site  !A
 year  !I
 pos   !I
 block  *
 row *
 plant !I
 clone !I 
 budset

Ekebo_Savar_budset.txt  !SKIP 1 !MAXIT 20 !DOPART 1 !DDF 2

!PART 1
!obtain data summary for sites
TABULATE budset ~ site !STATS
budset ~ mu site site.block site.year
         !r clone

#!PIN !DEFINE
#P Additive 1*2.0
#P Phenotypic 1+2
#H h2 3 4 
