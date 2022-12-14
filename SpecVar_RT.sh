#/bin/bash

# Relevant contexts identification by SpecVar, updated at Dec 14th, 2022 by Zhanying Feng, please contact zyfeng@amss.ac.cn to report bugs

PYTHON_HOME="/home/users/zyfeng/Software/anaconda2/bin/"
LDSC_HOME="/home/users/zyfeng/Software/ldsc/"
GWAS=$1
echo "Step 1: Estimating "${GWAS}"'s relevance to 77 human contexts"
$PYTHON_HOME/python $LDSC_HOME/munge_sumstats.py --sumstats ./Input/${GWAS}.txt --merge-alleles ./Data/w_hm3.snplist --out ./Input/${GWAS} --a1-inc
$PYTHON_HOME/python $LDSC_HOME/ldsc.py --h2 ./Input/${GWAS}.sumstats.gz --ref-ld-chr ./Data/SpecVar/SpecVar. --w-ld-chr ./Data/weights/weights. --overlap-annot --frqfile-chr ./Data/Frq/1000G.EAS.QC. --out ./Results/${GWAS}.SpecVar
cat ./Results/${GWAS}.SpecVar.results | sed '1,2d' | sed 's/_0//g' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' > ./Results/${GWAS}.SpecVar.RS.txt
rm -f ./Results/${GWAS}.SpecVar.log
sed s/GWAS/${GWAS}/g ./scr/GetRevTissue.py > GetRevTissue.py
$PYTHON_HOME/python GetRevTissue.py
rm -f GetRevTissue.py

echo "Step 2: Extracting "${GWAS}"'s SNP-associated regulatory network in relevant contexts"
for t in `cat ./Results/${GWAS}_SigTissue.txt`
do
	source AnnotSNP.sh $GWAS $t
done
