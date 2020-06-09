import numpy as np

f = open('Tissue_CRS_GWAS.txt')
CRS = f.readlines();f.close()
C = []
for i in range(len(CRS)):
	CRS[i] = CRS[i].split('\t')
	C.append(float(CRS[i][2]))
f = open('Tissue_RETG_GWAS_K.txt')
a = f.readlines();f.close()
K = [[],[]]
for i in range(len(a)):
	a[i] = a[i].split('\t')
	K[0].append(a[i][0]+'\t'+a[i][1])
	K[1].append(float(a[i][2]))
f = open('Tissue_RETG_GWAS_V.txt')
a = f.readlines();f.close()
V = [[],[],[]]
for i in range(len(a)):
	a[i] = a[i].split('\t')
	V[0].append(a[i][0])
	a[i][1] = a[i][1].split(' ')
	V[1].append(float(a[i][1][0]))
	V[2].append(a[i][1][1])

MaxC = max(C);MaxK = max(K[1]);MaxV = max(V[1])
f = open('./Input/GWAS.bed')
a = f.readlines();f.close()
SNPInfo = [[],[]]
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	SNPInfo[0].append(a[i][3])
	SNPInfo[1].append(a[i][4])

g = open('Tissue_RETG_GWAS_rs.txt','w')
for i in range(len(CRS)):
	indel1 = V[0].index(CRS[i][0])
	indel2 = K[0].index(CRS[i][0]+'\t'+CRS[i][1])
	CC = float(CRS[i][2]) / MaxC;KK = K[1][indel2] / MaxK;VV = V[1][indel1] / MaxV
	RS = (CC*KK*VV)**(1/3)
	g.write(CRS[i][0]+'\t'+CRS[i][1]+'\t'+str(CC)+'\t'+str(KK)+'\t'+str(VV)+'\t'+str(RS)+'\t')
	if 'NA' not in V[2][indel1]:
		TempV2 = V[2][indel1].strip('\n').split(';;')
		for j in range(len(TempV2)):
			TempV2[j] = TempV2[j].split(':')
			indel3 = SNPInfo[0].index(TempV2[j][0])
			TempV2[j] = TempV2[j][0]+'('+SNPInfo[1][indel3]+'):'+TempV2[j][1]
		TempV2 = ';;'.join(TempV2)
		g.write(TempV2+'\n')

	else:
		g.write(V[2][indel1])
g.close()
