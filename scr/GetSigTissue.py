import numpy as np

f = open('./Results/GWAS.PECA_specific.sorted.results')
a = f.readlines();f.close()
score = []
for j in range(len(a)):
	a[j] = a[j].split('\t')
	a[j][1] = float(a[j][1])
	score.append(a[j][1])
me = np.mean(score); sd = np.std(score)
cut = 2.325 * sd + me
SigTissue = []
for j in range(len(a)):
	if a[j][1] > cut:
		SigTissue.append(a[j][0])
if len(SigTissue) > 0:
	g = open('./Results/GWAS_SigTissue.txt','w')
	for j in range(len(SigTissue)):
		g.write(SigTissue[j]+'\n')
	g.close()
