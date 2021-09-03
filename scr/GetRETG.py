f = open('../Data/Networks/Tissue_network.txt')
a = f.readlines();f.close()
del a[0]
g = open('Tissue_RETG.txt','w')
RETG = []
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	if a[i][3] != '':
		a[i][3] = a[i][3].split(';')
		for j in range(len(a[i][3])):
			RETG.append(a[i][3][j].replace('_','\t')+'\t'+a[i][3][j]+'\t'+a[i][1]+'\n')
RETG = list(set(RETG))
RETG.sort()
for i in range(len(RETG)):
	g.write(RETG[i])
g.close()
