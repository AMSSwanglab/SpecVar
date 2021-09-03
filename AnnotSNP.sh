Trait=$1;t=$2
mkdir $Trait;cd $Trait
mkdir Results

echo "Annotating "${Trait}" in "${t}
\cp ../scr/* .
sed -i s/Tissue/${t}/g *py; sed -i s/GWAS/${Trait}/g *py;

python3 GetRETG.py
bedtools intersect -wa -wb -a ../Data/SpecificPeak/${t}.bed -b ${t}_RETG.txt | awk -v OFS='\t' '{print $5,$6,$7,$8,$9}' > ${t}_specific_RETG.bed
cat ../Input/${Trait}.bed | awk '$5<0.05' > ${Trait}.bed
bedtools slop -i ${Trait}.bed -g ../Data/hg19.sizes -b 50000 > a;cat ${Trait}.bed | awk '{print $2}' > b;paste a b > ${Trait}.bed;rm -f a b
bedtools intersect -wa -wb -a ${t}_specific_RETG.bed -b ./${Trait}.bed | awk -v OFS='\t' '{print $4"--"$5,$9"--"$10"--"$11"--"$12}' | sort -k1 > ${t}_specific_RETG_${Trait}.txt
source Merge.sh ${t}_specific_RETG_${Trait}.txt
cat ${t}_specific_RETG_${Trait}.txt | awk -F'--' '{print $1}' | sort | uniq | awk -F'_' -v OFS='\t' '{print $1,$2,$3,$0}' | sortBed > ${t}_specific_RE_${Trait}.txt
findMotifsGenome.pl ${t}_specific_RE_${Trait}.txt hg19 ./. -size given -find ../Data/all_motif_rmdup -preparsedDir ../Data/Homer/ -p 8 > ${t}_specific_${Trait}_MotifTarget.bed
cat ${t}_specific_${Trait}_MotifTarget.bed | awk 'NR>1' | cut -f 1,4,6 > ${t}_specific_${Trait}_MotifTarget.txt
rm -f ${t}_specific_${Trait}_MotifTarget.bed motifFindingParameters.txt

for RE in `cat ${t}_specific_RE_${Trait}.txt |awk '{print $4}'`
do
	cat ../Data/Networks/${t}_network.txt | grep ${RE} | awk -v OFS='\t' '{print $1,$2,$3}' | sort | uniq > .${t}_${Trait}_RE_TFTG.txt
	cat ../Data/openness/${t}_openness.bed | grep ${RE} > .${t}_${Trait}_RE_Opn.txt
	cat ../Data/openness/${t}_peak_gene_100k_corr.bed | grep ${RE} | awk -v OFS='\t' '{print $2,$4}' > .${t}_${Trait}_RE_Corr.txt
	cat ./${t}_specific_${Trait}_MotifTarget.txt | grep ${RE} | awk -v OFS='\t' '{print $2,$3}' > .${t}_${Trait}_RE_Motif.txt
	python3 GetCRS.py
	cat .${t}_${Trait}_RE_CRS.txt >> ${t}_specific_RETG_${Trait}_C.txt
	rm -f .${t}_${Trait}_RE_*
done
source Merge.sh ${t}_specific_RETG_${Trait}_C.txt
python3 MergeCRS.py
python3 GetK.py
python3 GetAS.py

mv Results/${Trait}_${t}_SubNetwork.txt ../Results/
cd ..
rm -rf $Trait
