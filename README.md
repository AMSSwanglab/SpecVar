# SpecVar
SpecVar is a convenient tool for estimating heritability enrichment, identifying relevant tissues and relevance correlation for human complex traits.

## Installation

1. Run the following commands for installation:<br>
```bash
  wget https://github.com/AMSSwanglab/SpecVar/archive/master.zip
    
  unzip master.zip
    
  cd SpecVar-master
```
2. Download the necessary files for SpecVar into **SpecVar-master** at: [https://drive.google.com/file/d/1mm_koaVDlYwrd3IRj76_EvC-iBWhvZC8/view?usp=sharing](https://drive.google.com/file/d/1mm_koaVDlYwrd3IRj76_EvC-iBWhvZC8/view?usp=sharing) and run the fowllowing command:<br>
```bash
  tar -zxvf Prior.tar.gz
```
3.  Install LDSC at: [https://github.com/bulik/ldsc](https://github.com/bulik/ldsc)<br>
Install Homer at: [http://homer.ucsd.edu/homer/download.html](http://homer.ucsd.edu/homer/download.html)<br>

4.  Edit the PATHON_HOME and LDSC_HOME path in SpecVar_GC.sh to your personal LDSC home to finish the installation

## Run relevant contexts idenfication mode of SpecVar: <br>
Taking GWAS of Educational Attainment (EA) and Cognitive Performance (CP) as example.<br>
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
2. After preparing the input files as above, run the following command to estimate genetic correlation:
```bash
  bash SpecVar_GC.sh EA CP
```
3. Output files in **Results** fold:<br>
**EA.SpecVar.RS.txt**: the EA's RS scores to 77 human contexts; <br> 
**EA_SigTissue.txt**: the EA's relevant human contexts; <br> 
**CP.SpecVar.RS.txt**: the CP's RS scores to 77 human contexts; <br> 
**CP_SigTissue.txt**: the CP's relevant human contexts; <br> 
**EA_CP_GC.txt**: the genetic correlation and p-value of EA and CP; <br>
**EA_frontal_cortex_SubNetwork.txt**: the EA's SNP associated regulatory subnetwork in the most relevant common context "frontal cortex"; <br>
**CP_frontal_cortex_SubNetwork.txt**: the CP's SNP associated regulatory subnetwork in the most relevant common context "frontal cortex". <br>


## Run relevant contexts idenfication mode of SpecVar: <br>
Taking GWAS of Educational Attainment (EA) and Cognitive Performance (CP) as example.<br>
1. Input files: EA.txt, EA.bed and CP.txt, CP.bed with same formats above <br>
2. After preparing the input files as above, run the following command to estimate genetic correlation:
```bash
  bash SpecVar_GC.sh EA CP
```
3. Output files in **Results** fold:<br>
**EA.SpecVar.RS.txt**: the EA's RS scores to 77 human contexts; <br> 
**EA_SigTissue.txt**: the EA's relevant human contexts; <br> 
**CP.SpecVar.RS.txt**: the CP's RS scores to 77 human contexts; <br> 
**CP_SigTissue.txt**: the CP's relevant human contexts; <br> 
**EA_CP_GC.txt**: the genetic correlation and p-value of EA and CP; <br>
**EA_frontal_cortex_SubNetwork.txt**: the EA's SNP associated regulatory subnetwork in the most relevant common context "frontal cortex"; <br>
**CP_frontal_cortex_SubNetwork.txt**: the CP's SNP associated regulatory subnetwork in the most relevant common context "frontal cortex". <br>

## Requirements

  Python <br>
  Python package: numpy, sklearn, and scipy <br>
  LDSC <br>
  Homer <br>
  Better if number of processors is more than 8. <br>
  
## Citation
If you use SpecVar or SpecVar associated concepts, please cite

Zhanying, et al. Cellular context-specific regulation-based inference of genetic correlation across human complex traits. 2021.
