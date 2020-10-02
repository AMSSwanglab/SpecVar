# vPECA2
vPECA2 is a convenient tool for Pinpointing relevant tissue and interpreting genetic variants with regulatory networks.

## Installation

1.  Install LDSC at: [https://github.com/bulik/ldsc](https://github.com/bulik/ldsc)<br>
Install Homer at: [http://homer.ucsd.edu/homer/download.html](http://homer.ucsd.edu/homer/download.html)<br>

2.  Run the following commands to finish installation:<br>
    <br>
    tar -jxvf Files.tar.bz2 <br>
    <br>
    wget https://github.com/AMSSwanglab/vPECA2/archive/master.zip <br>
    <br>
    unzip master.zip<br>
    <br>
    cd vPECA2-master

3.  Download the necessary files for vPECA2 into **vPECA2-master** at: <br>

    (1) [https://pan.baidu.com/s/1ydaXIA8naAUSWAwhEjZZIw](https://pan.baidu.com/s/1ydaXIA8naAUSWAwhEjZZIw) with extraction code: qq3o; Or <br>
    (2) [https://drive.google.com/file/d/1knVx0UjeMTGClgSCU8FxCFHO9Wg3_r-o/view?usp=sharing](https://drive.google.com/file/d/1knVx0UjeMTGClgSCU8FxCFHO9Wg3_r-o/view?usp=sharing) <br>

4.  Edit first line of vPECA2.sh to your personal LDSC home


## Run vPECA (Face GWAS as example)
```bash
bash vPECA2.sh Face
```

## Requirements

Python environment: python 3 <br>
Python package: numpy, sklearn, and scipy <br>
LDSC <br>
Homer <br>
