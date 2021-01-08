#/bin/bash

# SpecVar v4.0.1 updated Jan 8th 2021
# Step 1: Finding relevant tissues
# Step 2: Finding Associated regulatory netowrk

Trait=$1

echo "Step 1: Finding relevant tissues..."
python ~/software/ldsc/munge_sumstats.py --sumstats ./Input/${Trait}.txt --merge-alleles ./Data/w_hm3.snplist --out ./Input/${Trait} --a1-inc
python ~/software/ldsc/ldsc.py --h2 ./Input/${Trait}.sumstats.gz --ref-ld-chr ./Data/PECA_specific/SpecVar. --w-ld-chr ./Data/weights/weights. --overlap-annot --frqfile-chr ./Data/Frq/1000G.EAS.QC. --out ./Results/${Trait}.SpecVar
cat ./Results/${Trait}.SpecVar.results | sed '1,2d' | sed 's/_0//g' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' | sort -k2nr > ./Results/${Trait}.SpecVar.sorted.results
rm -f ./Results/${Trait}.SpecVar.log
sed s/GWAS/${Trait}/g ./scr/GetSigTissue.py >  GetSigTissue.py
python3 GetSigTissue.py
rm -f GetSigTissue.py

if [ `cat ./Results/${Trait}_SigTissue.txt | wc -l ` -gt 0 ]
then
	echo "Step 2: Finding Associated regulatory netowrk..."
	for t in `cat ./Results/${Trait}_SigTissue.txt`
	do
		echo "Annotating "${Trait}" in "${t}
		\cp ./scr/* .
		sed -i s/Tissue/${t}/g *py; sed -i s/GWAS/${Trait}/g *py;
	
		python3 GetRETG.py
		cat ./Input/${Trait}.bed | awk '$5<1/100000'  > ${Trait}_Sig.bed
		bedtools intersect -wa -wb -a ${t}_RETG.txt -b ${Trait}_Sig.bed | awk -v OFS='\t' '{print $4"--"$5,$9"--"$10"--"$11"--"$12"--"$13}' | sort -k1 > ${t}_RETG_${Trait}.txt
		if [ `cat ${t}_RETG_${Trait}.txt | wc -l ` -gt 0 ]
		then
			echo "Outputing network associated with significant SNPs......"
			source Merge.sh ${t}_RETG_${Trait}.txt	
			cat ${t}_RETG_${Trait}.txt | awk -F'--' '{print $1}' | sort | uniq | awk -F'_' -v OFS='\t' '{print $1,$2,$3,$0}' | sortBed > ${t}_RE_${Trait}.txt
			findMotifsGenome.pl ${t}_RE_${Trait}.txt hg19 ./. -size given -find ./Data/all_motif_rmdup -preparsedDir ./Data/Homer/ > ${t}_${Trait}_MotifTarget.bed
			cat ${t}_${Trait}_MotifTarget.bed | awk 'NR>1' | cut -f 1,4,6 > ${t}_${Trait}_MotifTarget.txt
			rm -f ${t}_${Trait}_MotifTarget.bed motifFindingParameters.txt

			for RE in `cat ${t}_RE_${Trait}.txt | awk '{print $4}'`
			do
				cat ./Data/Networks/${t}_network.txt | grep ${RE} | awk -v OFS='\t' '{print $1,$2,$3}' | sort | uniq > .${t}_${Trait}_RE_TFTG.txt
				cat ./Data/openness/${t}_openness.bed | grep ${RE} > .${t}_${Trait}_RE_Opn.txt
				cat ./Data/openness/${t}_peak_gene_100k_corr.bed | grep ${RE} | awk -v OFS='\t' '{print $2,$4}' > .${t}_${Trait}_RE_Corr.txt
				cat ./${t}_${Trait}_MotifTarget.txt | grep ${RE} | awk -v OFS='\t' '{print $2,$3}' > .${t}_${Trait}_RE_Motif.txt
				python3 GetCRS.py
				cat .${t}_${Trait}_RE_CRS.txt >> ${t}_RETG_${Trait}_C.txt
				rm -f .${t}_${Trait}_RE_*
			done
			source Merge.sh ${t}_RETG_${Trait}_C.txt
			python3 MergeCRS_Sig.py
			python3 GetK_Sig.py
			python3 GetRS_Sig.py
			rm -f ${Trait}_Sig.bed
		fi
		echo "Outputing specific network associated with SNPs......"
		bedtools intersect -wa -wb -a ./Data/PECA_specific_Peak/${t}.bed -b ${t}_RETG.txt | awk -v OFS='\t' '{print $5,$6,$7,$8,$9}' > ${t}_specific_RETG.bed
		bedtools intersect -wa -wb -a ${t}_specific_RETG.bed -b ./Input/${Trait}.bed | awk -v OFS='\t' '{print $4"--"$5,$9"--"$10"--"$11"--"$12"--"$13}' | sort -k1 > ${t}_specific_RETG_${Trait}.txt
		source Merge.sh ${t}_specific_RETG_${Trait}.txt
		cat ${t}_specific_RETG_${Trait}.txt | awk -F'--' '{print $1}' | sort | uniq | awk -F'_' -v OFS='\t' '{print $1,$2,$3,$0}' | sortBed > ${t}_specific_RE_${Trait}.txt
		findMotifsGenome.pl ${t}_specific_RE_${Trait}.txt hg19 ./. -size given -find ./Data/all_motif_rmdup -preparsedDir ./Data/Homer/ > ${t}_specific_${Trait}_MotifTarget.bed
		cat ${t}_specific_${Trait}_MotifTarget.bed | awk 'NR>1' | cut -f 1,4,6 > ${t}_specific_${Trait}_MotifTarget.txt
		rm -f ${t}_specific_${Trait}_MotifTarget.bed motifFindingParameters.txt

		for RE in `cat ${t}_specific_RE_${Trait}.txt |awk '{print $4}'`
		do
			cat ./Data/Networks/${t}_network.txt | grep ${RE} | awk -v OFS='\t' '{print $1,$2,$3}' | sort | uniq > .${t}_${Trait}_RE_TFTG.txt
			cat ./Data/openness/${t}_openness.bed | grep ${RE} > .${t}_${Trait}_RE_Opn.txt
			cat ./Data/openness/${t}_peak_gene_100k_corr.bed | grep ${RE} | awk -v OFS='\t' '{print $2,$4}' > .${t}_${Trait}_RE_Corr.txt
			cat ./${t}_specific_${Trait}_MotifTarget.txt | grep ${RE} | awk -v OFS='\t' '{print $2,$3}' > .${t}_${Trait}_RE_Motif.txt
			python3 GetCRS.py
			cat .${t}_${Trait}_RE_CRS.txt >> ${t}_specific_RETG_${Trait}_C.txt
			rm -f .${t}_${Trait}_RE_*
		done
		source Merge.sh ${t}_specific_RETG_${Trait}_C.txt
		python3 MergeCRS_Spe.py
		python3 GetK_Spe.py
		python3 GetRS_Spe.py
		rm -f ${t}*_RE* *py Merge.sh *MotifTarget*
	done
echo ${Trait}": Step 2 of SpecVar Done"
fi
