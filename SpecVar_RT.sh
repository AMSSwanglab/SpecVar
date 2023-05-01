#/bin/bash

# Relevant contexts identification version of SpecVar, updated at Dec 14th, 2022 by Zhanying Feng, please contact zyfeng@amss.ac.cn to report bugs

PYTHON_HOME="/home/users/zyfeng/Software/anaconda2/bin/"
PYTHON3_HOME="/home/users/zyfeng/Software/anaconda3/bin/"
LDSC_HOME="/home/users/zyfeng/Software/ldsc/"
Bedtools_HOME="/home/users/zyfeng/Software/bin/"
HOMER_HOME="/home/users/zyfeng/Software/bin/"
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
mkdir $GWAS;cd $GWAS
mkdir Results
for t in `cat ../Results/${GWAS}_SigTissue.txt | head -n 1`
do
	echo "Annotating "${GWAS}" in "${t}
	
	sed s/Tissue/${t}/g ../scr/GetRETG.py | sed s/GWAS/${GWAS}/g > GetRETG.py;
	$PYTHON3_HOME/python GetRETG.py;rm -f GetRETG.py
	$Bedtools_HOME/bedtools intersect -wa -wb -a ../Data/SpecificPeak/${t}.bed -b ${t}_RETG.txt | awk -v OFS='\t' '{print $5,$6,$7,$8,$9}' > ${t}_specific_RETG.bed
	cat ../Input/${GWAS}.bed | awk '$5<0.05' > ${GWAS}.bed
	$Bedtools_HOME/bedtools slop -i ${GWAS}.bed -g ../Data/hg19.sizes -b 50000 > a;cat ${GWAS}.bed | awk '{print $2}' > b;paste a b > ${GWAS}.bed;rm -f a b
	$Bedtools_HOME/bedtools intersect -wa -wb -a ${t}_specific_RETG.bed -b ./${GWAS}.bed | awk -v OFS='\t' '{print $4"--"$5,$9"--"$10"--"$11"--"$12}' | sort -k1 > ${t}_specific_RETG_${GWAS}.txt
	source ../scr/Merge.sh ${t}_specific_RETG_${GWAS}.txt
	cat ${t}_specific_RETG_${GWAS}.txt | awk -F'--' '{print $1}' | sort | uniq | awk -F'_' -v OFS='\t' '{print $1,$2,$3,$0}' | sortBed > ${t}_specific_RE_${GWAS}.txt
	$HOMER_HOME/findMotifsGenome.pl ${t}_specific_RE_${GWAS}.txt hg19 ./. -size given -find ../Data/all_motif_rmdup -preparsedDir ../Data/Homer/ -p 8 > ${t}_specific_${GWAS}_MotifTarget.bed
	cat ${t}_specific_${GWAS}_MotifTarget.bed | awk 'NR>1' | cut -f 1,4,6 > ${t}_specific_${GWAS}_MotifTarget.txt
	rm -f ${t}_specific_${GWAS}_MotifTarget.bed motifFindingParameters.txt

	sed s/Tissue/${t}/g ../scr/GetCRS.py | sed s/GWAS/${GWAS}/g > GetCRS.py
	for RE in `cat ${t}_specific_RE_${GWAS}.txt |awk '{print $4}'`
	do
		cat ../Data/Networks/${t}_network.txt | grep ${RE} | awk -v OFS='\t' '{print $1,$2,$3}' | sort | uniq > .${t}_${GWAS}_RE_TFTG.txt
		cat ../Data/openness/${t}_openness.bed | grep ${RE} > .${t}_${GWAS}_RE_Opn.txt
		cat ../Data/openness/${t}_peak_gene_100k_corr.bed | grep ${RE} | awk -v OFS='\t' '{print $2,$4}' > .${t}_${GWAS}_RE_Corr.txt
		cat ./${t}_specific_${GWAS}_MotifTarget.txt | grep ${RE} | awk -v OFS='\t' '{print $2,$3}' > .${t}_${GWAS}_RE_Motif.txt
		$PYTHON3_HOME/python GetCRS.py
		cat .${t}_${GWAS}_RE_CRS.txt >> ${t}_specific_RETG_${GWAS}_C.txt
		rm -f .${t}_${GWAS}_RE_*
	done
	rm -f GetCRS.py

	source ../scr/Merge.sh ${t}_specific_RETG_${GWAS}_C.txt
	sed s/Tissue/${t}/g ../scr/MergeCRS.py | sed s/GWAS/${GWAS}/g > MergeCRS.py;$PYTHON3_HOME/python MergeCRS.py;rm -f MergeCRS.py
	sed s/Tissue/${t}/g ../scr/GetK.py | sed s/GWAS/${GWAS}/g > GetK.py;$PYTHON3_HOME/python GetK.py;rm -f GetK.py
	sed s/Tissue/${t}/g ../scr/GetAS.py | sed s/GWAS/${GWAS}/g > GetAS.py;$PYTHON3_HOME/python GetAS.py;rm -f GetAS.py

	mv ./Results/${GWAS}_${t}_SubNetwork.txt ../Results/
done
cd ..
rm -rf $GWAS
