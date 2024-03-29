# SpecVar
SpecVar is a convenient tool for estimating heritability enrichment, identifying relevant tissues and relevance correlation for human complex traits.

## Installation

1. Run the following commands for installation:<br>
```bash
wget https://github.com/AMSSwanglab/SpecVar/archive/master.zip

unzip master.zip

cd SpecVar-master

wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1NNK_R_1DQo8wxkIoXLE2ntM4GilqmP60' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1NNK_R_1DQo8wxkIoXLE2ntM4GilqmP60" -O SpecVar_Data.tar.gz && rm -rf /tmp/cookies.txt

tar -zxvf SpecVar_Data.tar.gz

rm -f SpecVar_Data.tar.gz
```
Note: if you can't use wget command to download files in google drive, you can manually download SpecVar_Data.tar.gz at [https://drive.google.com/file/d/1NNK_R_1DQo8wxkIoXLE2ntM4GilqmP60/view?usp=share_link](https://drive.google.com/file/d/1NNK_R_1DQo8wxkIoXLE2ntM4GilqmP60/view?usp=share_link)<br>

2.  Install LDSC at: [https://github.com/bulik/ldsc](https://github.com/bulik/ldsc)<br>
Install Homer at: [http://homer.ucsd.edu/homer/download.html](http://homer.ucsd.edu/homer/download.html)<br>
Install bedtools at: [https://bedtools.readthedocs.io/en/latest/content/installation.html](https://bedtools.readthedocs.io/en/latest/content/installation.html)<br>

3.  Edit the **PATHON_HOME (python2), PATHON3_HOME (python3), LDSC_HOME, Bedtools_HOME, and HOMER_HOME** path in SpecVar_RT.sh, SpecVar_GC.sh, SpecVar_FE.sh to your personal software location to finish the installation

## Run relevant contexts idenfication mode of SpecVar: <br>
Taking GWAS of Educational Attainment (EA) as example.<br>
1. Input files: EA.txt, EA.bed<br>
```bash
  head EA.txt
  MarkerName      A1      A2      EAF     p       N
  rs13090388      C       T       0.6905  4.29e-54        1070751
  rs7630869       C       T       0.6922  4.61e-54        1070751
  rs7623659       T       C       0.3095  4.75e-54        1070751
  rs11922013      G       C       0.6905  5.94e-54        1070751
  rs9859556       G       T       0.6905  6.03e-54        1070751
  rs6779524       C       T       0.6905  6.30e-54        1070751
  rs9871380       A       G       0.3095  6.68e-54        1070751
  rs9878943       G       A       0.6905  6.68e-54        1070751
  rs9814873       G       A       0.3095  6.78e-54        1070751
  
  head EA.bed
  chr3    49391082        49391082        rs13090388      4.29e-54        ###
  chr3    49522543        49522543        rs7630869       4.61e-54        ###
  chr3    49414791        49414791        rs7623659       4.75e-54        ###
  chr3    49458355        49458355        rs11922013      5.94e-54        ###
  chr3    49455986        49455986        rs9859556       6.03e-54        ###
  chr3    49450449        49450449        rs6779524       6.30e-54        120.853
  chr3    49438221        49438221        rs9871380       6.68e-54        ###
  chr3    49434654        49434654        rs9878943       6.68e-54        ###
  chr3    49454112        49454112        rs9814873       6.78e-54        120.853
  chr3    49453834        49453834        rs6997  6.88e-54        120.853
```
2. After preparing the input files as above, run the following command to identify relevant contexts:
```bash
  bash SpecVar_RT.sh EA
```
3. Output files in **Results** fold:<br>
**EA.SpecVar.RS.txt**: the EA's RS scores to 77 human contexts; <br> 
**EA_SigTissue.txt**: the EA's relevant human contexts; <br> 
**EA_frontal_cortex_SubNetwork.txt**: the EA's SNP associated regulatory subnetwork in the most relevant context "frontal cortex"; <br>

## Run relevance correlation estimation mode of SpecVar: <br>
Taking GWAS of Educational Attainment (EA) and Cognitive Performance (CP) as example.<br>
1. Input files: EA.txt, EA.bed and CP.txt, CP.bed with same formats above <br>
2. After preparing the input files as above, run the following command to estimate phenotypic correlation:
```bash
  bash SpecVar_GC.sh EA CP
```
3. Output files in **./Results/** fold:<br>
**EA.SpecVar.RS.txt**: the EA's RS scores to 77 human contexts; <br> 
**EA_SigTissue.txt**: the EA's relevant human contexts; <br> 
**CP.SpecVar.RS.txt**: the CP's RS scores to 77 human contexts; <br> 
**CP_SigTissue.txt**: the CP's relevant human contexts; <br> 
**EA_CP_GC.txt**: the genetic correlation and p-value of EA and CP; <br>
**EA_frontal_cortex_SubNetwork.txt**: the EA's SNP associated regulatory subnetwork in the most relevant common context "frontal cortex"; <br>
**CP_frontal_cortex_SubNetwork.txt**: the CP's SNP associated regulatory subnetwork in the most relevant common context "frontal cortex". <br>

## Run fold enrichment mode of SpecVar: <br>
This program takes a set of genomic loci (SNPs or all kinds of QTLs) as input.Taking GWAS of COVID19 as example:
```bash
head COVID19.bed
chr12   103014757       103014757       12:103014757:C:A        1.652e-10
chr12   113352159       113352159       12:113352159:T:C        2.739e-08
chr12   113357193       113357193       12:113357193:G:A        3.097e-09
chr12   113357209       113357209       12:113357209:G:A        5.386e-09
chr12   113357442       113357442       12:113357442:G:A        5.575e-09
chr12   113359132       113359132       12:113359132:A:G        4.568e-08
chr12   113360737       113360737       12:113360737:G:T        2.933e-08
chr12   113361117       113361117       12:113361117:C:T        4.526e-09
chr12   113361158       113361158       12:113361158:C:T        4.889e-08
chr12   113361174       113361174       12:113361174:T:C        4.861e-08
```
Note: the last column contains the P-values of SNPs, which is required. If you don't have P-values for your loci, this can be replaced by any biologically meaningful scores or blank characters.<br>

Then run the following command to do FE analysis of SpecVar.
```bash
bash SpecVar_FE.sh COVID19
```
The Output files in **./Results/COVID19/** fold:<br>
**COVID19_FE.txt**: The COVID19's FE scores to 77 human contexts; <br> 
**COVID19_Primary_monocytes_from_peripheral_blood_network.txt**: The COVID19's associated regulatory network in "Primary_monocytes_from_peripheral_blood". Each row means one SNP with P-value is inside of RE, which will potentially affect TF's binding and TG's expression with regulatory strength TRS. <br>

## Requirements

  Python2 and Python3 <br>
  Python package: numpy, sklearn, and scipy <br>
  LDSC <br>
  Bedtools <br>
  Homer <br>
  Better if number of processors is more than 8. <br>
  
## Citation
If you use SpecVar or SpecVar associated concepts, please cite

Zhanying Feng, et al. [Heritability enrichment in context-specific regulatory networks improves phenotype-relevant tissue identification](https://elifesciences.org/articles/82535). eLife 2022.
