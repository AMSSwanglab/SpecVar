import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def GetPvalue1(X,x):
	nor = GaussianMixture(n_components=2).fit(X)
	mean = nor.means_;var = nor.covariances_;weight = nor.weights_
	norm1 = norm(loc=mean[0],scale=np.sqrt(var[0]));norm2 = norm(loc=mean[1],scale=np.sqrt(var[1]))
	p = 1 - ((norm1.cdf(x) * weight[0] + norm2.cdf(x) * weight[1])[0][0])
	return p

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

MaxC = max(max(Result1[3]),max(Result2[3]));MaxK = max(max(Result1[4]),max(Result2[4]))
AllRE1 = [];AllC1 = [];AllK1 = [];AllRS1 = [];AllSign1 = [];AllTF1 = []
AllRE2 = [];AllC2 = [];AllK2 = [];AllRS2 = [];AllSign2 = [];AllTF2 = []
for i in range(len(Result1[0])):
	CC = Result1[3][i]/MaxC;KK = Result1[4][i]/MaxK;RS = (CC*KK)**(1/2)
	if Result1[0][i] in Result2[0]:
		AllRE1.append(Result1[0][i]+'\t'+Result1[1][i]+'\t'+Result1[2][i])
		AllC1.append(CC);AllK1.append(KK);AllRS1.append(RS);AllSign1.append('1\t1');AllTF1.append(Result1[5][i])
	else:
		AllRE1.append(Result1[0][i]+'\t'+Result1[1][i]+'\t'+Result1[2][i])
		AllC1.append(CC);AllK1.append(KK);AllRS1.append(RS);AllSign1.append('1\t0');AllTF1.append(Result1[5][i])
for i in range(len(Result2[0])):
	CC = Result2[3][i]/MaxC;KK = Result2[4][i]/MaxK;RS = (CC*KK)**(1/2)
	if Result2[0][i] in Result1[0]:
		AllRE2.append(Result2[0][i]+'\t'+Result2[1][i]+'\t'+Result2[2][i])
		AllC2.append(CC);AllK2.append(KK);AllRS2.append(RS);AllSign2.append('1\t1');AllTF2.append(Result2[5][i])
	if Result2[0][i] not in Result1[0]:
		AllRE2.append(Result2[0][i]+'\t'+Result2[1][i]+'\t'+Result2[2][i])
		AllC2.append(CC);AllK2.append(KK);AllRS2.append(RS);AllSign2.append('0\t1');AllTF2.append(Result2[5][i])

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
mean2 = np.mean(AllRS2);sd2 = np.std(AllRS2)
def GetPvalue2(x):
	x = (x-mean2)/sd2
	return 1-norm.cdf(x)

g = open('./Results/GWAS_Tissue_Report.txt','w')
g.write('RE\tTG\tSNP\tRS\tP-value\tSigSNP\tSigRE\tTF\n')
for i in range(len(Result1[0])):
	pvalue=GetPvalue1(AllRS1[i])
	g.write(AllRE1[i]+'\t'+str(AllRS1[i])+'\t'+str(pvalue)+'\t'+AllSign1[i]+'\t'+AllTF1[i])
for i in range(len(Result2[0])):
	if Result2[0][i] not in Result1[0]:
		pvalue=GetPvalue2(AllRS2[i])
		g.write(AllRE2[i]+'\t'+str(AllRS2[i])+'\t'+str(pvalue)+'\t'+AllSign2[i]+'\t'+AllTF2[i])
g.close()

f = open('./Data/SpecificGene/TissueSpeGene.txt')
SpeTF = f.readlines();f.close()
SpeTF = [SpeTF[i].strip('\n') for i in range(len(SpeTF))]

f = open('./Results/GWAS_Tissue_Report.txt')
a = f.readlines();f.close()
del a[0]
g = open('Tissue_RETG_GWAS_Report_Clean.txt','w')
g.write('RE\tTG\tSNP\tRS\tP-Value\tSigSNP\tSpeRE\tTF\n')
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	a[i][7] = a[i][7].split(',')
	tmpTF = []
	for j in range(len(a[i][7])):
		if a[i][7][j] in SpeTF:
			tmpTF.append(a[i][7][j])
	if len(tmpTF) > 0:
		tmpTF = ','.join(tmpTF)
		g.write(a[i][0]+'\t'+a[i][1]+'\t'+a[i][2]+'\t'+a[i][3]+'\t'+a[i][4]+'\t'+a[i][5]+'\t'+a[i][6]+'\t'+tmpTF+'\n')
g.close()
