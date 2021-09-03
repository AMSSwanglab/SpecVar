import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

f = open('Tissue_specific_RETG_GWAS_C.txt')
CRS = f.readlines();f.close()

f = open('Tissue_specific_RETG_GWAS_K.txt')
a = f.readlines();f.close()
K = [[],[],[],[],[]]
for i in range(len(a)):
	a[i] = a[i].split('\t')
	K[0].append(a[i][0])
	K[1].append(a[i][2])
	K[2].append(a[i][3])
	K[3].append(a[i][4])
	K[4].append(float(a[i][5]))

Summary = [[],[],[],[],[],[],[],[]]
for i in range(len(CRS)):
	CRS[i] = CRS[i].split('\t')
	indel = K[0].index(CRS[i][0])
	Summary[0].append(CRS[i][0]);Summary[1].append(CRS[i][1]);
	Summary[2].append(K[1][indel]);Summary[3].append(K[2][indel]);Summary[4].append(K[3][indel]);
	Summary[5].append(float(CRS[i][2]));Summary[6].append(K[4][indel]);
	Summary[7].append(CRS[i][3])

MaxC = max(Summary[5]);MaxK = max(Summary[6])
AllRE = [];AllC = [];AllK = [];AllRS = [];AllTF = []
for i in range(len(Summary[0])):
	CC = Summary[5][i]/MaxC;KK = Summary[6][i]/MaxK;RS = (CC*KK)**(1/2)
	AllRE.append(Summary[0][i]+'\t'+Summary[1][i]+'\t'+Summary[2][i]+'\t'+Summary[3][i]+'\t'+Summary[4][i])
	AllC.append(CC);AllK.append(KK);AllRS.append(RS);AllTF.append(Summary[7][i])

mean = np.mean(AllRS);sd = np.std(AllRS)
def GetPvalue2(x):
	x = (x-mean)/sd
	return 1-norm.cdf(x)

g = open('./Results/GWAS_Tissue_SubNetwork.txt','w')
g.write('RE\tTG\tSNP\tD0\tSig\tRS\tP-value\tTF\n')
for i in range(len(Summary[0])):
	indel = AllRS.index(max(AllRS))
	pvalue=GetPvalue2(AllRS[indel])
	g.write(AllRE[indel]+'\t'+str(AllRS[indel])+'\t'+str(pvalue)+'\t'+AllTF[indel])
	AllRS[indel] = -100
g.close()
