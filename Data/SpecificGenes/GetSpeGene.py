import numpy as np

f = open('../../../V4/Comparison/LDSC/Genes/Exp.txt')
a = f.readlines();f.close()
head = a[0].strip('\n').split('\t');del head[0]
del a[0]
gene = []
for i in range(len(a)):
	a[i] = a[i].split('\t')
	gene.append(a[i][0]);del a[i][0]
	for j in range(len(a[i])):
		a[i][j] = float(a[i][j])
	mea = np.mean(a[i]);sd = np.std(a[i])
	for j in range(len(a[i])):
		if sd > 0:
			a[i][j] = (a[i][j]-mea)/sd
		else:
			a[i][j] = 0.0
for i in range(len(head)):
	g = open(head[i]+'SpeGene.txt','w')
	for j in range(len(a)):
		if a[j][i] >= 1.5:
			g.write(gene[j]+'\n')
	g.close()
