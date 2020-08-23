import numpy as np

f = open('RE_TFTG.txt')
a = f.readlines();f.close()
TFTG = [[],[],[]]
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	TFTG[0].append(a[i][0])
	TFTG[1].append(a[i][1])
	TFTG[2].append(float(a[i][2]))
f = open('RE_Corr.txt')
a = f.readlines();f.close()
peak_corr = [[],[]]
for i in range(len(a)):
	a[i] = a[i].split('\t')
	peak_corr[0].append(a[i][0])
	peak_corr[1].append(float(a[i][1]))
f = open('./Data/TFMotifMatch.txt')
a = f.readlines();f.close()
TFMotif = [[],[]]
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	TFMotif[0].append(a[i][0]);del a[i][0]
	TFMotif[1].append(a[i])
f = open('RE_Motif.txt')
a = f.readlines();f.close()
REMotif = [[],[]]
for i in range(len(a)):
	a[i] = a[i].split('\t')
	REMotif[0].append(a[i][0])
	REMotif[1].append(float(a[i][1]))
f = open('RE_Opn.txt')
opn = f.readline();f.close()
opn = opn.split('\t');opn[1] = float(opn[1])

TG = list(set(TFTG[1]))
g = open('RE_CRS.txt','w')
for i in range(len(TG)):
	CRS1 = 0
	CRS_TF = []
	for j in range(len(TFTG[0])):
		if TFTG[1][j] == TG[i]:
			Bik = []
			Motif = TFMotif[1][TFMotif[0].index(TFTG[0][j])]
			for k in range(len(Motif)):
				if Motif[k] in REMotif[0]:
					Bik.append(REMotif[1][REMotif[0].index(Motif[k])])
			if len(Bik) > 0:
				CRS_TF.append(TFTG[0][j])
				Bik = np.max(Bik)
				CRS1 += Bik*TFTG[2][j]
			else:
				CRS1 += 0
	CRS2 = np.max([np.abs(peak_corr[1][j]) for j in range(len(peak_corr[0])) if peak_corr[0][j]==TG[i]])
	CRS = CRS1*CRS2*opn[1]
	CRS_TF = '__'.join(list(set(CRS_TF)))
	g.write(opn[0]+'\t'+TG[i]+'--'+str(CRS)+'--'+CRS_TF+'\n')
g.close()
