import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

f = open('Tissue_specific_RETG_GWAS_C.txt')
CRS = f.readlines();f.close()

f = open('Tissue_specific_RETG_GWAS_K.txt')
a = f.readlines();f.close()
K = [[],[],[]]
for i in range(len(a)):
	a[i] = a[i].split('\t')
	K[0].append(a[i][0])
	K[1].append(a[i][2])
	K[2].append(float(a[i][3]))

Result2 = [[],[],[],[],[],[]]
for i in range(len(CRS)):
	CRS[i] = CRS[i].split('\t')
	indel = K[0].index(CRS[i][0])
	Result2[0].append(CRS[i][0]);Result2[1].append(CRS[i][1]);Result2[2].append(K[1][indel]);Result2[3].append(float(CRS[i][2]));
	Result2[4].append(K[2][indel]);Result2[5].append(CRS[i][3])

MaxC = max(Result2[3]);MaxK = max(Result2[4])
AllRE2 = [];AllC2 = [];AllK2 = [];AllRS2 = [];AllTF2 = []
for i in range(len(Result2[0])):
	CC = Result2[3][i]/MaxC;KK = Result2[4][i]/MaxK;RS = (CC*KK)**(1/2)
	AllRE2.append(Result2[0][i]+'\t'+Result2[1][i]+'\t'+Result2[2][i])
	AllC2.append(CC);AllK2.append(KK);AllRS2.append(RS);AllTF2.append(Result2[5][i])

mean2 = np.mean(AllRS2);sd2 = np.std(AllRS2)
def GetPvalue2(x):
	x = (x-mean2)/sd2
	return 1-norm.cdf(x)

g = open('./Results/GWAS_Tissue_Report.txt','w')
g.write('RE\tTG\tSNP\tRS\tP-value\tTF\n')
for i in range(len(Result2[0])):
	indel = AllRS2.index(max(AllRS2))
	pvalue=GetPvalue2(AllRS2[indel])
	g.write(AllRE2[indel]+'\t'+str(AllRS2[indel])+'\t'+str(pvalue)+'\t'+AllTF2[indel])
	AllRS2[indel] = -100
g.close()

f = open('./Data/SpecificGene/TissueSpeGene.txt')
SpeTF = f.readlines();f.close()
SpeTF = [SpeTF[i].strip('\n') for i in range(len(SpeTF))]

f = open('./Results/GWAS_Tissue_Report.txt')
a = f.readlines();f.close()
del a[0]
g = open('./Results/GWAS_Tissue_Report_Clean.txt','w')
g.write('RE\tTG\tSNP\tRS\tP-Value\tTF\n')
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	a[i][5] = a[i][5].split(',')
	tmpTF = []
	for j in range(len(a[i][5])):
		if a[i][5][j] in SpeTF:
			tmpTF.append(a[i][5][j])
	if len(tmpTF) > 0:
		tmpTF = ','.join(tmpTF)
		g.write(a[i][0]+'\t'+a[i][1]+'\t'+a[i][2]+'\t'+a[i][3]+'\t'+a[i][4]+'\t'+tmpTF+'\n')
g.close()
