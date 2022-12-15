#/bin/bash

# Phenotypic Correlation version of SpecVar 1.0, updated at Dec 14th, 2022 by Zhanying Feng, please contact zyfeng@amss.ac.cn to report bugs

PYTHON_HOME="~/software/anaconda2/bin/"
LDSC_HOME='~/software/ldsc'
#Step 1: Estimating traits' relevance to 77 human contexts

GWAS1=$1;GWAS2=$2
echo "Estimating "${GWAS1}"'s relevance to 77 human contexts"
$PYTHON_HOME/python $LDSC_HOME/munge_sumstats.py --sumstats ./Input/${GWAS1}.txt --merge-alleles ./Data/w_hm3.snplist --out ./Input/${GWAS1} --a1-inc
$PYTHON_HOME/python $LDSC_HOME/ldsc.py --h2 ./Input/${GWAS1}.sumstats.gz --ref-ld-chr ./Data/SpecVar/SpecVar. --w-ld-chr ./Data/weights/weights. --overlap-annot --frqfile-chr ./Data/Frq/1000G.EAS.QC. --out ./Results/${GWAS1}.SpecVar
cat ./Results/${GWAS1}.SpecVar.results | sed '1,2d' | sed 's/_0//g' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' > ./Results/${GWAS1}.SpecVar.RS.txt
rm -f ./Results/${GWAS1}.SpecVar.log
sed s/GWAS/${GWAS1}/g ./scr/GetRevTissue.py > GetRevTissue.py
$PYTHON_HOME/python GetRevTissue.py
rm -f GetRevTissue.py

echo "Estimating "${GWAS2}"'s relevance to 77 human contexts"
$PYTHON_HOME/python $LDSC_HOME/munge_sumstats.py --sumstats ./Input/${GWAS2}.txt --merge-alleles ./Data/w_hm3.snplist --out ./Input/${GWAS2} --a1-inc
$PYTHON_HOME/python $LDSC_HOME/ldsc.py --h2 ./Input/${GWAS2}.sumstats.gz --ref-ld-chr ./Data/SpecVar/SpecVar. --w-ld-chr ./Data/weights/weights. --overlap-annot --frqfile-chr ./Data/Frq/1000G.EAS.QC. --out ./Results/${GWAS2}.SpecVar
cat ./Results/${GWAS2}.SpecVar.results | sed '1,2d' | sed 's/_0//g' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' > ./Results/${GWAS2}.SpecVar.RS.txt
rm -f ./Results/${GWAS2}.SpecVar.log
sed s/GWAS/${GWAS2}/g ./scr/GetRevTissue.py > GetRevTissue.py
$PYTHON_HOME/python GetRevTissue.py
rm -f GetRevTissue.py

cat ./scr/GetGC.py | sed s/GWAS1/${GWAS1}/g | sed s/GWAS2/${GWAS2}/g > GetGC.py
$PYTHON_HOME/python GetGC.py
rm -f GetGC.py

cat ./scr/GetMRC.py | sed s/GWAS1/${GWAS1}/g | sed s/GWAS2/${GWAS2}/g > GetMRC.py
MRC=`$PYTHON_HOME/python GetMRC.py`
rm -f GetMRC.py
echo "The Most Relevant Common Context is: "${MRC}

source AnnotSNP.sh $GWAS1 $MRC
source AnnotSNP.sh $GWAS2 $MRC
echo "SpecVar for "${GWAS1}" and "${GWAS2}" done!"
