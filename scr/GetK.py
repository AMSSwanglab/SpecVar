import numpy as np

f = open('Tissue_RETG_GWAS.txt')
a = f.readlines();f.close()
g = open('Tissue_RETG_GWAS_K.txt','w')
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	a[i][1] = a[i][1].split(',')
	K = 0
	for j in range(len(a[i][1])):
		a[i][1][j] = a[i][1][j].split('--')
		if a[i][1][j][4] == '###':
			K += (-np.log(float(a[i][1][j][1])))
		else:
			K += (-np.log(float(a[i][1][j][1]))/float(a[i][1][j][4]))
	g.write(a[i][0].replace('--','\t')+'\t'+str(K)+'\n')
g.close()
