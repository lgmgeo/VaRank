#!/bin/bash


mkdir -p $VARANK/Annotations_Exomiser/2007
cd $VARANK/Annotations_Exomiser/2007

# download
echo -e "\t...downloading 2007_hg19.tar.gz"
curl -C - -LO https://www.lbgi.fr/~geoffroy/Annotations/2007_hg19.tar.gz >& download-1.log

echo -e "\t...donwnloading 2007_phenotype.zip"
curl -C - -LO https://data.monarchinitiative.org/exomiser/data/2007_phenotype.zip >& download-2.log

# search for error
grep -i "error" download-1.log
grep -i "error" download-2.log

# uncompress
tar -xf 2007_hg19.tar.gz >& uncompress-1.log
unzip 2007_phenotype.zip >& uncompress-2.log

# search for error
grep -i "error" uncompress-1.log
grep -i "error" uncompress-2.log

# clean
rm -rf 2007_phenotype.zip
rm -rf 2007_hg19.tar.gz
rm -f 2007_phenotype.sha256
rm -f download-1.log
rm -f download-2.log

# finished
echo -e "\t...done ($VARANK/Annotations_Exomiser/2007/*)"

