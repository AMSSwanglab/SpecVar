import numpy as np

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
g = open('Tissue_specific_RETG_GWAS.txt','w')
for i in range(len(RE)):
	g.write(RE[i]+'--'+','.join(TG[i])+'\t'+SNP[i])
g.close()
f = open('Tissue_specific_RETG_GWAS.txt')
a = f.readlines();f.close()
g = open('Tissue_specific_RETG_GWAS_K.txt','w')
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	a[i][1] = a[i][1].split(',')
	K = 0;KK = []
	for j in range(len(a[i][1])):
		a[i][1][j] = a[i][1][j].split('--')
		if a[i][1][j][4] == '###':
			K += (-np.log(float(a[i][1][j][1])))
			KK.append(a[i][1][j][0])
		else:
			K += (-np.log(float(a[i][1][j][1]))/float(a[i][1][j][4]))
			KK.append(a[i][1][j][0])
	K = K / len(a[i][1])
	g.write(a[i][0].replace('--','\t')+'\t'+','.join(KK)+'\t'+str(K)+'\n')
g.close()
