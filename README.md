# vPECA2
vPECA2 is a convenient tool for Pinpointing relevant tissue and interpreting genetic variants with regulatory networks.

## Installation

1.  Install LDSC at: [https://github.com/bulik/ldsc](https://github.com/bulik/ldsc)<br>
Install Homer at: [http://homer.ucsd.edu/homer/download.html](http://homer.ucsd.edu/homer/download.html)<br>

2.  Run the following commands to start installation:<br>
```bash
  tar -jxvf Files.tar.bz2
    
  wget https://github.com/AMSSwanglab/vPECA2/archive/master.zip
    
  unzip master.zip
    
  cd vPECA2-master
```
3.  Download the necessary files for vPECA2 into **vPECA2-master** at: <br>
    (1) [https://pan.baidu.com/s/1WHmyg06Ob6XXCLfmXJl-IA](https://pan.baidu.com/s/1WHmyg06Ob6XXCLfmXJl-IA) with extraction code: 11ih ; Or <br>
    (2) [https://drive.google.com/file/d/17Rrysp64sS0tum4WB0ONcxthmDGBEKbw/view?usp=sharing](https://drive.google.com/file/d/17Rrysp64sS0tum4WB0ONcxthmDGBEKbw/view?usp=sharing) <br>

4.  Edit first line of vPECA2.sh to your personal LDSC home to finish the installation


## Run vPECA (Face GWAS as example)
```bash
  bash vPECA2.sh Face
```

## Requirements

Python environment: python 3 <br>
Python package: numpy, sklearn, and scipy <br>
LDSC <br>
Homer <br>
