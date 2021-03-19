# VaRank: An integrated tool for SNV/indel annotation and ranking 
https://lbgi.fr/VaRank/

## QUICK INSTALLATION

1. The sources can be cloned to any directory:
```
cd /path/to/install/
git clone git@github.com:lgmgeo/VaRank.git
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
- $ALAMUT: Alamut Batch installation directory
- $SNPEFF: SnpEff and SnpSift installation directory

## TEST with the Alamut annotation engine

1. Change to the repo directory, look at the configfile and run the test
```
cd $VARANK/Tests/Alamut/

cat configfile

$VARANK/bin/VaRank -vcfDir "." -SamOut "Sample1 Sample2" 
```
2. Examine the outputs (*tsv)

