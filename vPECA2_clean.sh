python ~/software/ldsc/munge_sumstats.py --sumstats ./Input/${1}.txt --merge-alleles ~/software/ldsc/${1}/w_hm3.snplist --out ./Input/${1} --a1-inc

python ~/software/ldsc/ldsc.py --h2 ./Input/${1}.sumstats.gz --ref-ld-chr ./Data/PECA_specific/86Tissue. --w-ld-chr ~/software/ldsc/weights/weights. --overlap-annot --frqfile-chr ./Data/Frq/1000G.EAS.QC. --out ./Result/${1}.PECA_specific

cat ./Result/${1}.PECA_specific.results | sed '1,2d' | sed 's/_0//g' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' | sort -k2nr > ./Result/${1}.PECA_specific.sorted.results

echo "The Top 2 tissue associated with "${1}" is:"
head -n 2 ./Result/${1}.PECA_specific.sorted.results
head -n 2 ./Result/${1}.PECA_specific.sorted.results | awk '{print $1}' | sed s/_0//g > ${1}_Top2_Tissue.txt

for t in `head -n 1 ./Result/${1}_Top2_Tissue.txt`
do
	echo "Annotating "${1}" in "${t}
	\cp ./scr/* .
	sed -i s/Tissue/${t}/g *py; sed -i s/GWAS/${1}/g *py;
	
	echo "Step 1: Annotating "${1}" significant SNPs..."

	python3 GetRETG.py
	bedtools intersect -wa -wb -a ${t}_RETG.txt -b ./Input/${1}.bed | awk -v OFS='\t' '{print $4"--"$5,$9"--"$10"--"$11"--"$12"--"$13}' | sort -k1 > ${t}_RETG_${1}.txt
	source Merge.sh ${t}_RETG_${1}.txt
		
	cat ${t}_RETG_${1}.txt | awk -F'--' '{print $1}' | sort | uniq > ${t}_RE_${1}.txt
	for RE in `cat ${t}_RE_${1}.txt`
	do
		cat ./Data/Networks/${t}_network.txt | grep ${RE} | awk -v OFS='\t' '{print $1,$2,$3}' | sort | uniq > RE_TFTG.txt
		cat ./Data/openness/${t}_openness.bed | grep ${RE} > RE_Opn.txt
		cat ./Data/openness/${t}_peak_gene_100k_corr.bed | grep ${RE} | awk -v OFS='\t' '{print $2,$4}' > RE_Corr.txt
		cat ./Data/MotifTarget/${t}_MotifTarget.txt | grep ${RE} | awk -v OFS='\t' '{print $2,$3}' > RE_Motif.txt
		python3 GetCRS.py
		cat RE_CRS.txt >> ${t}_RETG_${1}_C.txt
		rm -f RE_*
	done
	source Merge.sh ${t}_RETG_${1}_C.txt
	python3 MergeCRS.py
	python3 GetK.py
	
	echo "Step 2: Annotating "${t}" PECA2 specific REs..."
	
	bedtools intersect -wa -wb -a ./Data/PECA_specific_Peak/${t}.bed -b ${t}_RETG.txt | awk -v OFS='\t' '{print $5,$6,$7,$8,$9}' > ${t}_specific_RETG.bed
	bedtools intersect -wa -wb -a ${t}_specific_RETG.bed -b ./Input/${1}_full.bed | awk -v OFS='\t' '{print $4"--"$5,$9"--"$10"--"$11"--"$12"--"$13}' | sort -k1 > ${t}_specific_RETG_${1}.txt
	source Merge.sh ${t}_specific_RETG_${1}.txt
	cat ${t}_specific_RETG_${1}.txt | awk -F'--' '{print $1}' | sort | uniq > ${t}_specific_RE_${1}.txt
	for RE in `cat ${t}_specific_RE_${1}.txt`
	do
		cat ./Data/Networks/${t}_network.txt | grep ${RE} | awk -v OFS='\t' '{print $1,$2,$3}' | sort | uniq > RE_TFTG.txt
		cat ./Data/openness/${t}_openness.bed | grep ${RE} > RE_Opn.txt
		cat ./Data/openness/${t}_peak_gene_100k_corr.bed | grep ${RE} | awk -v OFS='\t' '{print $2,$4}' > RE_Corr.txt
		cat ./Data/MotifTarget/${t}_MotifTarget.txt | grep ${RE} | awk -v OFS='\t' '{print $2,$3}' > RE_Motif.txt
		python3 GetCRS.py
		cat RE_CRS.txt >> ${t}_specific_RETG_${1}_C.txt
		rm -f RE_*
	done
	source Merge.sh ${t}_specific_RETG_${1}_C.txt
	python3 MergeCRS_full.py
	python3 GetK_full.py

	python3 GetRS.py
	echo -e -n 'RE\tTG\tC\tK\tRS\tSigSNP\tSpeRE\tSignificant\tTF\n' > ${1}_${t}_Report_Clean_Sorted.txt
	sed '1d' ${t}_RETG_${1}_Report_Clean.txt | sort -k6nr,6 -k7nr,7 -k5nr,5 >> ${1}_${t}_Report_Clean_Sorted.txt
	rm -f ${t}*_RE* *py Merge.sh
done
