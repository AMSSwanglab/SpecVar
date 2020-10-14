LDSC_HOME='/home/fengzhanying//software/ldsc'
Homer_HOME='/home/fengzhanying/software/PECA/PECA-master-v3'

python ${LDSC_HOME}/munge_sumstats.py --sumstats ./Input/${1}.txt --merge-alleles ./Data/w_hm3.snplist --out ./Input/${1} --a1-inc

echo "Step 1: Finding relevant tissues..."
python ${LDSC_HOME}/ldsc.py --h2 ./Input/${1}.sumstats.gz --ref-ld-chr ./Data/PECA_specific/86Tissue. --w-ld-chr ./Data/weights/weights. --overlap-annot --frqfile-chr ./Data/Frq/1000G.EAS.QC. --out ./Results/${1}.PECA_specific

cat ./Results/${1}.PECA_specific.results | sed '1,2d' | sed 's/_0//g' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' | sort -k2nr > ./Results/${1}.PECA_specific.sorted.results
rm -f ./Results/${1}.PECA_specific.results
rm -f ./Results/${1}.PECA_specific.log

sed s/GWAS/${1}/g ./scr/GetSigTissue.py >  GetSigTissue.py
python3  GetSigTissue.py

echo "Step 2: Finding Associated regulatory netowrk..."
for t in `cat ./Results/${1}_SigTissue.txt`
do
	echo "Annotating "${1}" in "${t}
	\cp ./scr/* .
	sed -i s/Tissue/${t}/g *py; sed -i s/GWAS/${1}/g *py;
	
	python3 GetRETG.py
	bedtools intersect -wa -wb -a ${t}_RETG.txt -b ./Input/${1}.bed | awk -v OFS='\t' '{print $4"--"$5,$9"--"$10"--"$11"--"$12"--"$13}' | sort -k1 > ${t}_RETG_${1}.txt
	if [ `cat ${t}_RETG_${1}.txt | wc -l ` -gt 0 ]
	then
		source Merge.sh ${t}_RETG_${1}.txt	
		cat ${t}_RETG_${1}.txt | awk -F'--' '{print $1}' | sort | uniq | awk -F'_' -v OFS='\t' '{print $1,$2,$3,$0}' | sortBed > ${t}_RE_${1}.txt
		findMotifsGenome.pl ${t}_RE_${1}.txt hg19 ./. -size given -find ./Data/all_motif_rmdup -preparsedDir ${Homer_HOME}/Homer/ > ${t}_${1}_MotifTarget.bed
		cat ${t}_${1}_MotifTarget.bed | awk 'NR>1' | cut -f 1,4,6 > ${t}_${1}_MotifTarget.txt
		rm -f ${t}_${1}_MotifTarget.bed motifFindingParameters.txt

		for RE in `cat ${t}_RE_${1}.txt | awk '{print $4}'`
		do
			cat ./Data/Networks/${t}_network.txt | grep ${RE} | awk -v OFS='\t' '{print $1,$2,$3}' | sort | uniq > RE_TFTG.txt
			cat ./Data/openness/${t}_openness.bed | grep ${RE} > RE_Opn.txt
			cat ./Data/openness/${t}_peak_gene_100k_corr.bed | grep ${RE} | awk -v OFS='\t' '{print $2,$4}' > RE_Corr.txt
			cat ./${t}_${1}_MotifTarget.txt | grep ${RE} | awk -v OFS='\t' '{print $2,$3}' > RE_Motif.txt
			python3 GetCRS.py
			cat RE_CRS.txt >> ${t}_RETG_${1}_C.txt
			rm -f RE_*
		done
		source Merge.sh ${t}_RETG_${1}_C.txt
		python3 MergeCRS.py
		python3 GetK.py
	
		bedtools intersect -wa -wb -a ./Data/PECA_specific_Peak/${t}.bed -b ${t}_RETG.txt | awk -v OFS='\t' '{print $5,$6,$7,$8,$9}' > ${t}_specific_RETG.bed
		bedtools intersect -wa -wb -a ${t}_specific_RETG.bed -b ./Input/${1}_full.bed | awk -v OFS='\t' '{print $4"--"$5,$9"--"$10"--"$11"--"$12"--"$13}' | sort -k1 > ${t}_specific_RETG_${1}.txt
		source Merge.sh ${t}_specific_RETG_${1}.txt
		cat ${t}_specific_RETG_${1}.txt | awk -F'--' '{print $1}' | sort | uniq | awk -F'_' -v OFS='\t' '{print $1,$2,$3,$0}' | sortBed > ${t}_specific_RE_${1}.txt
		findMotifsGenome.pl ${t}_specific_RE_${1}.txt hg19 ./. -size given -find ./Data/all_motif_rmdup -preparsedDir ${Homer_HOME}/Homer/ > ${t}_specific_${1}_MotifTarget.bed
		cat ${t}_specific_${1}_MotifTarget.bed | awk 'NR>1' | cut -f 1,4,6 > ${t}_specific_${1}_MotifTarget.txt
		rm -f ${t}_specific_${1}_MotifTarget.bed motifFindingParameters.txt

		for RE in `cat ${t}_specific_RE_${1}.txt |awk '{print $4}'`
		do
			cat ./Data/Networks/${t}_network.txt | grep ${RE} | awk -v OFS='\t' '{print $1,$2,$3}' | sort | uniq > RE_TFTG.txt
			cat ./Data/openness/${t}_openness.bed | grep ${RE} > RE_Opn.txt
			cat ./Data/openness/${t}_peak_gene_100k_corr.bed | grep ${RE} | awk -v OFS='\t' '{print $2,$4}' > RE_Corr.txt
			cat ./${t}_specific_${1}_MotifTarget.txt | grep ${RE} | awk -v OFS='\t' '{print $2,$3}' > RE_Motif.txt
			python3 GetCRS.py
			cat RE_CRS.txt >> ${t}_specific_RETG_${1}_C.txt
			rm -f RE_*
		done
		source Merge.sh ${t}_specific_RETG_${1}_C.txt
		python3 MergeCRS_full.py
		python3 GetK_full.py

		python3 GetRS.py
		echo -e -n 'RE\tTG\tSNP\tRS\tP-value\tSigSNP\tSpeRE\tTF\n' > ./Results/${1}_${t}_Report_Clean_Sorted.txt
		sed '1d' ${t}_RETG_${1}_Report_Clean.txt | sort -k6nr,6 -k4nr,4 >> ./Results/${1}_${t}_Report_Clean_Sorted.txt
		rm -f ${t}*_RE* *py Merge.sh *MotifTarget*
	else
		bedtools intersect -wa -wb -a ./Data/PECA_specific_Peak/${t}.bed -b ${t}_RETG.txt | awk -v OFS='\t' '{print $5,$6,$7,$8,$9}' > ${t}_specific_RETG.bed
		bedtools intersect -wa -wb -a ${t}_specific_RETG.bed -b ./Input/${1}_full.bed | awk -v OFS='\t' '{print $4"--"$5,$9"--"$10"--"$11"--"$12"--"$13}' | sort -k1 > ${t}_specific_RETG_${1}.txt
		source Merge.sh ${t}_specific_RETG_${1}.txt
		cat ${t}_specific_RETG_${1}.txt | awk -F'--' '{print $1}' | sort | uniq | awk -F'_' -v OFS='\t' '{print $1,$2,$3,$0}' | sortBed > ${t}_specific_RE_${1}.txt
		findMotifsGenome.pl ${t}_specific_RE_${1}.txt hg19 ./. -size given -find ./Data/all_motif_rmdup -preparsedDir ${Homer_HOME}/Homer/ > ${t}_specific_${1}_MotifTarget.bed
		cat ${t}_specific_${1}_MotifTarget.bed | awk 'NR>1' | cut -f 1,4,6 > ${t}_specific_${1}_MotifTarget.txt
		rm -f ${t}_specific_${1}_MotifTarget.bed motifFindingParameters.txt

		for RE in `cat ${t}_specific_RE_${1}.txt |awk '{print $4}'`
		do
			cat ./Data/Networks/${t}_network.txt | grep ${RE} | awk -v OFS='\t' '{print $1,$2,$3}' | sort | uniq > RE_TFTG.txt
			cat ./Data/openness/${t}_openness.bed | grep ${RE} > RE_Opn.txt
			cat ./Data/openness/${t}_peak_gene_100k_corr.bed | grep ${RE} | awk -v OFS='\t' '{print $2,$4}' > RE_Corr.txt
			cat ./${t}_specific_${1}_MotifTarget.txt | grep ${RE} | awk -v OFS='\t' '{print $2,$3}' > RE_Motif.txt
			python3 GetCRS.py
			cat RE_CRS.txt >> ${t}_specific_RETG_${1}_C.txt
			rm -f RE_*
		done
		source Merge.sh ${t}_specific_RETG_${1}_C.txt
		python3 MergeCRS_full.py
		python3 GetK_full.py

		python3 GetRS2.py
		echo -e -n 'RE\tTG\tSNP\tRS\tP-value\tSigSNP\tSpeRE\tTF\n' > ./Results/${1}_${t}_Report_Clean_Sorted.txt
		sed '1d' ${t}_RETG_${1}_Report_Clean.txt | sort -k6nr,6 -k4nr,4 >> ./Results/${1}_${t}_Report_Clean_Sorted.txt
		rm -f ${t}*_RE* *py Merge.sh *MotifTarget*
	fi
done
echo ${1}" vPECA2 Done"
