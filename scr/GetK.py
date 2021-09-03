import numpy as np

d0 = 10000
f = open('Tissue_specific_RETG_GWAS.txt')
a = f.readlines();f.close()
RE = [];TG = [];SNP = []
for i in range(len(a)):
	a[i] = a[i].split('\t')
	a[i][0] = a[i][0].split('--')
	if a[i][0][0] not in RE:
		RE.append(a[i][0][0])
		TG.append([a[i][0][1]])
		SNP.append(a[i][1])
	else:
		indel = RE.index(a[i][0][0])
		TG[indel].append(a[i][0][1])
g = open('Tissue_specific_RETG_GWAS2.txt','w')
for i in range(len(RE)):
	g.write(RE[i]+'--'+','.join(TG[i])+'\t'+SNP[i])
g.close()
f = open('Tissue_specific_RETG_GWAS2.txt')
a = f.readlines();f.close()
g = open('Tissue_specific_RETG_GWAS_K.txt','w')
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	a[i][1] = a[i][1].split(',')
	K = 0;KK = [];K0 = [[],[]];KS = [[],[]]
	st = int(a[i][0].split('--')[0].split('_')[1]);ed = int(a[i][0].split('--')[0].split('_')[2])
	for j in range(len(a[i][1])):
		a[i][1][j] = a[i][1][j].split('--')
		if float(a[i][1][j][1]) == 0:
			a[i][1][j][1] = (0.1)**305
		distance = -1
		if (int(a[i][1][j][3])-st)*(int(a[i][1][j][3])-ed) < 0:
			distance = 0
			K0[0].append(a[i][1][j][0]);K0[1].append(float(a[i][1][j][1]))
		else:
			distance = min(abs(int(a[i][1][j][3])-st),abs(int(a[i][1][j][3])-ed))
		if float(a[i][1][j][1]) <= 5/100000000:
			KS[0].append(a[i][1][j][0]);KS[1].append(distance)
		
		if a[i][1][j][2] == '###':
			a[i][1][j][2] = 1.0
		K += (-np.log(float(a[i][1][j][1]))/float(a[i][1][j][2]))*np.e**(-distance/10000)
		KK.append(a[i][1][j][0])
	K = K / len(a[i][1])
	if len(K0[0]) == 0:
		K0 = "###"
	else:
		K0T = []
		for j in range(len(K0[0])):
			indel = K0[1].index(min(K0[1]))
			K0T.append(K0[0][indel]);K0[1][indel] = 100000000000000000000000
			
		K0 = ','.join(K0T)
	if len(KS[0]) == 0:
		KS = "###"
	else:
		KST = []
		for j in range(len(KS[1])):
			indel = KS[1].index(min(KS[1]))
			KST.append(KS[0][indel]);KS[1][indel] = 100000000000000000000000
		KS = ','.join(KST)
	g.write(a[i][0].replace('--','\t')+'\t'+','.join(KK)+'\t'+K0+'\t'+KS+'\t'+str(K)+'\n')
g.close()
