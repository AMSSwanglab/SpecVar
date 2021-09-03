# SpecVar
SpecVar is a convenient tool for estimating interpretable genetic correlation of human complex traits and annotating the SNPs with context specific regulatory networks.

## Installation

1. Run the following commands for installation:<br>
```bash
  tar -zxvf Prior.tar.gz
    
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

4.  Edit the LDSC_HOME path in SpecVar_GC.sh to your personal LDSC home to finish the installation


## Run SpecVar: <br>
Taking GWAS of Educational Attainment (EA) and Cognitive Performance (CP) as example.<br>
1. Input files: EA.txt, EA.bed and CP.txt, CP.bed <br>
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
```
```bash
  bash SpecVar_GC.sh EA CP
```
## Results (GWAS of Face as example)

**Face.SpecVar.results:** The LDSC standard output of partitioned heritability enrichment with 77 specific regulatory network. <br>
**Face.SpecVar.sorted.results:** Tissues ranked by RS score. <br>
**Face_SigTissue.txt:** The significantly relevant tissues. <br>
**Face_CNCC_Report.txt:** The associated CNCC-specific regulatory sub-network with Face GWAS. <br>
**Face_CNCC_SigSNP.txt:** The ssociated CNCC regulatory sub-network with significant SNPs of Face GWAS. <br>

## Requirements

  Python <br>
  Python package: numpy, sklearn, and scipy <br>
  LDSC <br>
  Homer <br>
  Better if number of processors is more than 8. <br>
  
## Citation
