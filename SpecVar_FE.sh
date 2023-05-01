#/bin/bash

# Fold enrichment version of SpecVar 1.0, updated at May 1st, 2023 by Zhanying Feng, please contact zyfeng@amss.ac.cn to report bugs

GWAS=$1
mkdir ${GWAS}_RE_Overlap
for i in `cat ./Data/RELength.txt | awk '{print $1}'`
do
	bedtools intersect -wa -wb -a ./Input/${GWAS}.bed -b ./Data/PECA_Peaks/${i}.bed > ${GWAS}_RE_Overlap/${GWAS}_${i}.txt
	num=`bedtools intersect -wa -wb -a ./Input/${GWAS}.bed -b ./Data/PECA_Peaks/${i}.bed | awk '{print $4}' | sort | uniq | wl`
	echo -e -n ${i}"\t"${num}"\n" >> ${GWAS}_RE_Overlap.txt
done
sed s/GWAS/${GWAS}/g ./scr/GetFE.py > GetFE.py
mkdir ./Results/${GWAS}
python3 GetFE.py;rm -rf ${GWAS}_RE_Overlap ${GWAS}_RE_Overlap.txt GetFE.py
