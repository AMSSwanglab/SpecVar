python ~/software/ldsc/munge_sumstats.py --sumstats ./Input/${1}.txt --merge-alleles ~/software/ldsc/${1}/w_hm3.snplist --out ./Input/${1} --a1-inc

python ~/software/ldsc/ldsc.py --h2 ./Input/${1}.sumstats.gz --ref-ld-chr ./Data/Peak_full/86Tissue. --w-ld-chr ~/software/ldsc/weights/weights. --overlap-annot --frqfile-chr ./Data/Frq/1000G.EAS.QC. --out ./Result/${1}.Peak_full
cat ./Result/${1}.Peak_full.results | sed '1,2d' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' | sort -k2nr > ./Result/${1}.Peak_full.sorted.results
python ~/software/ldsc/ldsc.py --h2 ./Input/${1}.sumstats.gz --ref-ld-chr ./Data/Peak_specific/86Tissue. --w-ld-chr ~/software/ldsc/weights/weights. --overlap-annot --frqfile-chr ./Data/Frq/1000G.EAS.QC. --out ./Result/${1}.Peak_specific
cat ./Result/${1}.Peak_specific.results | sed '1,2d' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' | sort -k2nr > ./Result/${1}.PECA_full.sorted.results
python ~/software/ldsc/ldsc.py --h2 ./Input/${1}.sumstats.gz --ref-ld-chr ./Data/PECA_full/86Tissue. --w-ld-chr ~/software/ldsc/weights/weights. --overlap-annot --frqfile-chr ./Data/Frq/1000G.EAS.QC. --out ./Result/${1}.PECA_full
cat ./Result/${1}.PECA_full.results | sed '1,2d' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' | sort -k2nr > ./Result/${1}.PECA_full.sorted.results
python ~/software/ldsc/ldsc.py --h2 ./Input/${1}.sumstats.gz --ref-ld-chr ./Data/PECA_specific/86Tissue. --w-ld-chr ~/software/ldsc/weights/weights. --overlap-annot --frqfile-chr ./Data/Frq/1000G.EAS.QC. --out ./Result/${1}.PECA_specific
cat ./Result/${1}.PECA_specific.results | sed '1,2d' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' | sort -k2nr > ./Result/${1}.PECA_specific.sorted.results
python ~/software/ldsc/ldsc.py --h2 ./Input/${1}.sumstats.gz --ref-ld-chr ./Data/Gene_specific/86Tissue. --w-ld-chr ~/software/ldsc/weights/weights. --overlap-annot --frqfile-chr ./Data/Frq/1000G.EAS.QC. --out ./Result/${1}.Gene_specific

cat ./Result/${1}.Gene_specific.results | sed '1,2d' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' | sort -k2nr > ./Result/${1}.Gene_specific.sorted.results
cat ./Result/${1}.Peak_full.results | sed '1,2d' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' | sort -k2nr > ./Result/${1}.Peak_full.sorted.results
cat ./Result/${1}.Peak_specific.results | sed '1,2d' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' | sort -k2nr > ./Result/${1}.Peak_specific.sorted.results
cat ./Result/${1}.PECA_full.results | sed '1,2d' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' | sort -k2nr > ./Result/${1}.PECA_full.sorted.results
cat ./Result/${1}.PECA_specific.results | sed '1,2d' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' | sort -k2nr > ./Result/${1}.PECA_specific.sorted.results
cat ./Result/${1}.Gene_specific.results | sed '1,2d' | grep -v HumanBrain | awk -v OFS='\t' '{print $1,-$5*log($7)/log(10)}' | sort -k2nr > ./Result/${1}.Gene_specific.sorted.results

echo "The Top 2 tissue associated with "${1}" is:"
head -n 2 ./Result/${1}.PECA_specific.sorted.results | awk '{print $1}' | sed s/_0//g
head -n 2 ./Result/${1}.PECA_specific.sorted.results | awk '{print $1}' | sed s/_0//g > ${1}_Top2_Tissue.txt

for t in `cat ${1}_Top2_Tissue.txt`
do
	\cp ./scr/* .
	sed -i s/Tissue/${t}/g *py; sed -i s/GWAS/${1}/g *py;
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
		cat RE_CRS.txt >> ${t}_CRS_${1}.txt
		rm -f RE_*
	done
	
	python3 GetK.py
	
	for RE in `cat ${t}_RETG_${1}.txt | awk -F'--' '{print $1}'`
	do
		echo $RE | tr '_' '\t' > a.bed
		bedtools intersect -wa -wb -a ./Input/${1}.bed -b a.bed > a_snp.txt
		chr=`echo $RE | awk -F'_' '{print $1}'`
		st=`echo $RE | awk -F'_' '{print $2-1}'`
		ed=`echo $RE | awk -F'_' '{print $3}'`
		echo -e -n ${chr}'\t'${st}'\t'${ed} > a.bed
		bedtools getfasta -fi /home/fengzhanying/data/GenomeIndexFile/hg19.fa -bed a.bed > ${RE}.fa
		findMotifs.pl ${RE}.fa fasta ./. -find ./Data/all_motif_rmdup > ${RE}_res.txt
		rm -f *tmp motifFindingParameters.txt 
		python3 GetSNPFasta.py
		for snp in `cat a_snp.txt | awk '{print $4}'`
		do
			findMotifs.pl ${RE}_v_${snp}.fa fasta ./. -find ./Data/all_motif_rmdup > ${RE}_v_${snp}_res.txt
			rm -f *tmp motifFindingParameters.txt
		done
		V=`python3 GetVar.py`
		echo -e -n ${RE}'\t'$V'\n' >> ${t}_RETG_${1}_V.txt
		rm -f a.bed a_snp.txt *res.txt *Var.txt *fa
	done
	
	python3 GetRS.py
	sort -k6nr ${t}_RETG_${1}_rs.txt > ${1}_${t}_rs_sorted.txt
	python3 FilterTF.py
	rm -f ${t}_RE* ${t}_CRS* *py Merge.sh
done
