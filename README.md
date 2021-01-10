# SpecVar
SpecVar is a convenient tool for pinpointing relevant tissue and interpreting genetic variants with specific regulatory networks.

## Installation

1.  Install LDSC at: [https://github.com/bulik/ldsc](https://github.com/bulik/ldsc)<br>
Install Homer at: [http://homer.ucsd.edu/homer/download.html](http://homer.ucsd.edu/homer/download.html)<br>

2.  Run the following commands to start installation:<br>
```bash
  tar -jxvf Files.tar.bz2
    
  wget https://github.com/AMSSwanglab/SpecVar/archive/master.zip
    
  unzip master.zip
    
  cd SpecVar-master
```
3.  Download the necessary files for SpecVar into **SpecVar-master** at: <br>
    (1) [https://pan.baidu.com/s/1WHmyg06Ob6XXCLfmXJl-IA](https://pan.baidu.com/s/1WHmyg06Ob6XXCLfmXJl-IA) with extraction code: 11ih ; Or <br>
    (2) [https://drive.google.com/file/d/17Rrysp64sS0tum4WB0ONcxthmDGBEKbw/view?usp=sharing](https://drive.google.com/file/d/17Rrysp64sS0tum4WB0ONcxthmDGBEKbw/view?usp=sharing) <br>

4.  Edit first line of SpecVar.sh to your personal LDSC home to finish the installation


## Run SpecVar (GWAS of Facail distances (Face) as example)
```bash
  bash SpecVar.sh Face
```
## Results (GWAS of Face as example)

**Face.SpecVar.results:** The LDSC standard output of partitioned heritability enrichment with 77 specific regulatory network. <br>
**Face.SpecVar.sorted.results:** Tissues ranked by RS score. <br>
**Face_SigTissue.txt:** The significantly relevant tissues. <br>
**Face_CNCC_Report.txt:** The associated CNCC-specific regulatory sub-network with Face GWAS. <br>
**Face_CNCC_SigSNP.txt:** The ssociated CNCC regulatory sub-network with significant SNPs of Face GWAS. <br>

## Requirements

  Python environment: python 3 <br>
  Python package: numpy, sklearn, and scipy <br>
  LDSC <br>
  Homer <br>
  
## Citation
