#/bin/bash

# Phenotypic Correlation version of SpecVar 1.0, updated at Dec 14th, 2022 by Zhanying Feng, please contact zyfeng@amss.ac.cn to report bugs

PYTHON_HOME="/home/users/zyfeng/Software/anaconda2/bin/"
PYTHON3_HOME="/home/users/zyfeng/Software/anaconda3/bin/"
LDSC_HOME="/home/users/zyfeng/Software/ldsc/"
Bedtools_HOME="/home/users/zyfeng/Software/bin/"
HOMER_HOME="/home/users/zyfeng/Software/bin/"
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
$PYTHON3_HOME/python GetRevTissue.py
rm -f GetRevTissue.py

cat ./scr/GetGC.py | sed s/GWAS1/${GWAS1}/g | sed s/GWAS2/${GWAS2}/g > GetGC.py
$PYTHON3_HOME/python GetGC.py
rm -f GetGC.py

cat ./scr/GetMRC.py | sed s/GWAS1/${GWAS1}/g | sed s/GWAS2/${GWAS2}/g > GetMRC.py
MRC=`$PYTHON3_HOME/python GetMRC.py`
rm -f GetMRC.py
echo "The Most Relevant Common Context is: "${MRC}

echo "Annotating "${GWAS1}" in "${MRC}
mkdir $GWAS1;cd $GWAS1;mkdir Results

sed s/Tissue/${MRC}/g ../scr/GetRETG.py | sed s/GWAS/${GWAS1}/g > GetRETG.py;
$PYTHON3_HOME/python GetRETG.py;rm -f GetRETG.py
$Bedtools_HOME/bedtools intersect -wa -wb -a ../Data/SpecificPeak/${MRC}.bed -b ${MRC}_RETG.txt | awk -v OFS='\t' '{print $5,$6,$7,$8,$9}' > ${MRC}_specific_RETG.bed
cat ../Input/${GWAS1}.bed | awk '$5<0.05' > ${GWAS1}.bed
$Bedtools_HOME/bedtools slop -i ${GWAS1}.bed -g ../Data/hg19.sizes -b 50000 > a;cat ${GWAS1}.bed | awk '{print $2}' > b;paste a b > ${GWAS1}.bed;rm -f a b
$Bedtools_HOME/bedtools intersect -wa -wb -a ${MRC}_specific_RETG.bed -b ./${GWAS1}.bed | awk -v OFS='\t' '{print $4"--"$5,$9"--"$10"--"$11"--"$12}' | sort -k1 > ${MRC}_specific_RETG_${GWAS1}.txt
source ../scr/Merge.sh ${MRC}_specific_RETG_${GWAS1}.txt
cat ${MRC}_specific_RETG_${GWAS1}.txt | awk -F'--' '{print $1}' | sort | uniq | awk -F'_' -v OFS='\t' '{print $1,$2,$3,$0}' | sortBed > ${MRC}_specific_RE_${GWAS1}.txt
$HOMER_HOME/findMotifsGenome.pl ${MRC}_specific_RE_${GWAS1}.txt hg19 ./. -size given -find ../Data/all_motif_rmdup -preparsedDir ../Data/Homer/ -p 8 > ${MRC}_specific_${GWAS1}_MotifTarget.bed
cat ${MRC}_specific_${GWAS1}_MotifTarget.bed | awk 'NR>1' | cut -f 1,4,6 > ${MRC}_specific_${GWAS1}_MotifTarget.txt
rm -f ${MRC}_specific_${GWAS1}_MotifTarget.bed motifFindingParameters.txt

sed s/Tissue/${MRC}/g ../scr/GetCRS.py | sed s/GWAS/${GWAS1}/g > GetCRS.py
for RE in `cat ${MRC}_specific_RE_${GWAS1}.txt |awk '{print $4}'`
do
	cat ../Data/Networks/${MRC}_network.txt | grep ${RE} | awk -v OFS='\t' '{print $1,$2,$3}' | sort | uniq > .${MRC}_${GWAS1}_RE_TFTG.txt
	cat ../Data/openness/${MRC}_openness.bed | grep ${RE} > .${MRC}_${GWAS1}_RE_Opn.txt
	cat ../Data/openness/${MRC}_peak_gene_100k_corr.bed | grep ${RE} | awk -v OFS='\t' '{print $2,$4}' > .${MRC}_${GWAS1}_RE_Corr.txt
	cat ./${MRC}_specific_${GWAS1}_MotifTarget.txt | grep ${RE} | awk -v OFS='\t' '{print $2,$3}' > .${MRC}_${GWAS1}_RE_Motif.txt
	$PYTHON3_HOME/python GetCRS.py
	cat .${MRC}_${GWAS1}_RE_CRS.txt >> ${MRC}_specific_RETG_${GWAS1}_C.txt
	rm -f .${MRC}_${GWAS1}_RE_*
done
rm -f GetCRS.py

source ../scr/Merge.sh ${MRC}_specific_RETG_${GWAS1}_C.txt
sed s/Tissue/${MRC}/g ../scr/MergeCRS.py | sed s/GWAS/${GWAS1}/g > MergeCRS.py;$PYTHON3_HOME/python MergeCRS.py;rm -f MergeCRS.py
sed s/Tissue/${MRC}/g ../scr/GetK.py | sed s/GWAS/${GWAS1}/g > GetK.py;$PYTHON3_HOME/python GetK.py;rm -f GetK.py
sed s/Tissue/${MRC}/g ../scr/GetAS.py | sed s/GWAS/${GWAS1}/g > GetAS.py;$PYTHON3_HOME/python GetAS.py;rm -f GetAS.py

mv ./Results/${GWAS1}_${MRC}_SubNetwork.txt ../Results/
cd ..
rm -rf $GWAS1

echo "Annotating "${GWAS2}" in "${MRC}
mkdir $GWAS2;cd $GWAS2;mkdir Results

sed s/Tissue/${MRC}/g ../scr/GetRETG.py | sed s/GWAS/${GWAS2}/g > GetRETG.py;
$PYTHON3_HOME/python GetRETG.py;rm -f GetRETG.py
$Bedtools_HOME/bedtools intersect -wa -wb -a ../Data/SpecificPeak/${MRC}.bed -b ${MRC}_RETG.txt | awk -v OFS='\t' '{print $5,$6,$7,$8,$9}' > ${MRC}_specific_RETG.bed
cat ../Input/${GWAS2}.bed | awk '$5<0.05' > ${GWAS2}.bed
$Bedtools_HOME/bedtools slop -i ${GWAS2}.bed -g ../Data/hg19.sizes -b 50000 > a;cat ${GWAS2}.bed | awk '{print $2}' > b;paste a b > ${GWAS2}.bed;rm -f a b
$Bedtools_HOME/bedtools intersect -wa -wb -a ${MRC}_specific_RETG.bed -b ./${GWAS2}.bed | awk -v OFS='\t' '{print $4"--"$5,$9"--"$10"--"$11"--"$12}' | sort -k1 > ${MRC}_specific_RETG_${GWAS2}.txt
source ../scr/Merge.sh ${MRC}_specific_RETG_${GWAS2}.txt
cat ${MRC}_specific_RETG_${GWAS2}.txt | awk -F'--' '{print $1}' | sort | uniq | awk -F'_' -v OFS='\t' '{print $1,$2,$3,$0}' | sortBed > ${MRC}_specific_RE_${GWAS2}.txt
$HOMER_HOME/findMotifsGenome.pl ${MRC}_specific_RE_${GWAS2}.txt hg19 ./. -size given -find ../Data/all_motif_rmdup -preparsedDir ../Data/Homer/ -p 8 > ${MRC}_specific_${GWAS2}_MotifTarget.bed
cat ${MRC}_specific_${GWAS2}_MotifTarget.bed | awk 'NR>1' | cut -f 1,4,6 > ${MRC}_specific_${GWAS2}_MotifTarget.txt
rm -f ${MRC}_specific_${GWAS2}_MotifTarget.bed motifFindingParameters.txt

sed s/Tissue/${MRC}/g ../scr/GetCRS.py | sed s/GWAS/${GWAS2}/g > GetCRS.py
for RE in `cat ${MRC}_specific_RE_${GWAS2}.txt |awk '{print $4}'`
do
	cat ../Data/Networks/${MRC}_network.txt | grep ${RE} | awk -v OFS='\t' '{print $1,$2,$3}' | sort | uniq > .${MRC}_${GWAS2}_RE_TFTG.txt
	cat ../Data/openness/${MRC}_openness.bed | grep ${RE} > .${MRC}_${GWAS2}_RE_Opn.txt
	cat ../Data/openness/${MRC}_peak_gene_100k_corr.bed | grep ${RE} | awk -v OFS='\t' '{print $2,$4}' > .${MRC}_${GWAS2}_RE_Corr.txt
	cat ./${MRC}_specific_${GWAS2}_MotifTarget.txt | grep ${RE} | awk -v OFS='\t' '{print $2,$3}' > .${MRC}_${GWAS2}_RE_Motif.txt
	$PYTHON3_HOME/python GetCRS.py
	cat .${MRC}_${GWAS2}_RE_CRS.txt >> ${MRC}_specific_RETG_${GWAS2}_C.txt
	rm -f .${MRC}_${GWAS2}_RE_*
done
rm -f GetCRS.py

source ../scr/Merge.sh ${MRC}_specific_RETG_${GWAS2}_C.txt
sed s/Tissue/${MRC}/g ../scr/MergeCRS.py | sed s/GWAS/${GWAS2}/g > MergeCRS.py;$PYTHON3_HOME/python MergeCRS.py;rm -f MergeCRS.py
sed s/Tissue/${MRC}/g ../scr/GetK.py | sed s/GWAS/${GWAS2}/g > GetK.py;$PYTHON3_HOME/python GetK.py;rm -f GetK.py
sed s/Tissue/${MRC}/g ../scr/GetAS.py | sed s/GWAS/${GWAS2}/g > GetAS.py;$PYTHON3_HOME/python GetAS.py;rm -f GetAS.py

mv ./Results/${GWAS2}_${MRC}_SubNetwork.txt ../Results/
cd ..
rm -rf $GWAS2
