import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

f = open('Tissue_RETG_GWAS_C.txt')
CRS = f.readlines();f.close()

f = open('Tissue_RETG_GWAS_K.txt')
a = f.readlines();f.close()
K = [[],[],[]]
for i in range(len(a)):
	a[i] = a[i].split('\t')
	K[0].append(a[i][0])
	K[1].append(a[i][2])
	K[2].append(float(a[i][3]))

Result1 = [[],[],[],[],[],[]]
for i in range(len(CRS)):
	CRS[i] = CRS[i].split('\t')
	indel = K[0].index(CRS[i][0])
	Result1[0].append(CRS[i][0]);Result1[1].append(CRS[i][1]);Result1[2].append(K[1][indel]);Result1[3].append(float(CRS[i][2]));
	Result1[4].append(K[2][indel]);Result1[5].append(CRS[i][3])

MaxC = max(Result1[3]);MaxK = max(Result1[4])
AllRE1 = [];AllC1 = [];AllK1 = [];AllRS1 = [];AllTF1 = []
for i in range(len(Result1[0])):
	CC = Result1[3][i]/MaxC;KK = Result1[4][i]/MaxK;RS = (CC*KK)**(1/2)
	AllRE1.append(Result1[0][i]+'\t'+Result1[1][i]+'\t'+Result1[2][i])
	AllC1.append(CC);AllK1.append(KK);AllRS1.append(RS);AllTF1.append(Result1[5][i])

if len(AllRS1) > 20:
	nor = GaussianMixture(n_components=2).fit(np.array(AllRS1).reshape(len(AllRS1),1))
	mean = nor.means_;var = nor.covariances_;weight = nor.weights_
	norm1 = norm(loc=mean[0],scale=np.sqrt(var[0]));norm2 = norm(loc=mean[1],scale=np.sqrt(var[1]))
	def GetPvalue1(x):
		p = 1 - ((norm1.cdf(x) * weight[0] + norm2.cdf(x) * weight[1])[0][0])
		return p
else:
	def GetPvalue1(x):
		return 0

g = open('./Results/GWAS_Tissue_SigSNP.txt','w')
g.write('RE\tTG\tSNP\tRS\tP-value\tTF\n')
for i in range(len(Result1[0])):
	indel = AllRS1.index(max(AllRS1))
	pvalue=GetPvalue1(AllRS1[indel])
	g.write(AllRE1[indel]+'\t'+str(AllRS1[indel])+'\t'+str(pvalue)+'\t'+AllTF1[indel])
	AllRS1[indel] = -100
g.close()

f = open('./Data/SpecificGene/TissueSpeGene.txt')
SpeTF = f.readlines();f.close()
SpeTF = [SpeTF[i].strip('\n') for i in range(len(SpeTF))]

f = open('./Results/GWAS_Tissue_SigSNP.txt')
a = f.readlines();f.close()
del a[0]
g = open('./Results/GWAS_Tissue_SigSNP_Clean.txt','w')
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
