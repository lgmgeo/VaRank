# VaRank: An integrated tool for SNV/indel annotation and ranking 
https://lbgi.fr/VaRank/

## QUICK INSTALLATION

1. The sources can be cloned to any directory:
```
cd /path/to/install/
git clone https://github.com/lgmgeo/VaRank.git
```

2. Set the VARANK global environmental variable as the location of the git repo on your system. 

In csh:
```
setenv VARANK /path/to/install/VaRank/
```
In bash:
```
export VARANK=/path/to/install/VaRank
```
3. Set the ALAMUT/SNPEFF global environmental variable.
Depending on the selected annotation engine:
- $ALAMUT : Alamut Batch installation directory
- $SNPEFF : SnpEff and SnpSift installation directory

## TEST

1. Change to the repo directory, and run the test
```
cd /path/to/install/VaRank/bin/VaRank/SnpEffTests/
$VARANK/bin/VaRank -vcfdir "."
```
2. Examine the outputs (fam1_Sample1_*)

