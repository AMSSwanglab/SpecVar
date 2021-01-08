import numpy as np

f = open('Tissue_specific_RETG_GWAS_C.txt')
a = f.readlines();f.close()
g = open('Tissue_specific_RETG_GWAS_C.txt','w')
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	a[i][1]  = a[i][1].split(',')
	TG = [];CRS = [];TF = []
	for j in range(len(a[i][1])):
		a[i][1][j] = a[i][1][j].split('--')
		TG.append(a[i][1][j][0])
		CRS.append(float(a[i][1][j][1]))
		TF += a[i][1][j][2].split('__')
	g.write(a[i][0]+'\t'+','.join(TG)+'\t'+str(np.log2(np.mean(CRS)+1))+'\t'+','.join(list(set(TF)))+'\n')
g.close()
