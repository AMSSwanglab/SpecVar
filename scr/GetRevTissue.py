import numpy as np

f = open('./Results/GWAS.SpecVar.RS.txt')
a = f.readlines();f.close()
score = []
for j in range(len(a)):
	a[j] = a[j].split('\t')
	a[j][1] = float(a[j][1])
	score.append(a[j][1])
me = np.mean(score); sd = np.std(score)
cut = 1.65 * sd + me
SigTissue = [[],[]]
for j in range(len(a)):
	if a[j][1] > cut:
		SigTissue[0].append(a[j][0])
		SigTissue[1].append(a[j][1])
g = open('./Results/GWAS_SigTissue.txt','w')
if len(SigTissue) > 0:
	for j in range(len(SigTissue[0])):
		indel = SigTissue[1].index(max(SigTissue[1]))
		g.write(SigTissue[0][indel]+'\n')
		SigTissue[1][indel] = -10000
	g.close()
else:
	g.write(a[score.index(max(score))][0]+'\n')
