import numpy as np

f = open('./Data/MotifTFMatch.txt')
a = f.readlines();f.close()
MotifTF = [[],[]]
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	MotifTF[0].append(a[i][0])
	del a[i][0]
	MotifTF[1].append(','.join(a[i]))

f = open('a_snp.txt')
RE_SNP = f.readlines();f.close()
MaxV = [[],[],[]]
for i in range(len(RE_SNP)):
	RE_SNP[i] = RE_SNP[i].strip('\n').split('\t')
	RE = RE_SNP[i][8]+'_'+RE_SNP[i][9]+'_'+RE_SNP[i][10]
	f = open(RE+'_res.txt')
	a = f.readlines();f.close()
	O = [[],[],[]]
	del a[0]
	for j in range(len(a)):
		a[j] = a[j].split('\t')
		if a[j][3] in MotifTF[0]:
			O[0].append(a[j][1]+'\t'+a[j][3])
			O[1].append(float(a[j][5]))
			O[2].append(MotifTF[1][MotifTF[0].index(a[j][3])])
	f = open(RE+'_v_'+RE_SNP[i][3]+'_res.txt')
	a = f.readlines();f.close()
	N = [[],[],[]]
	del a[0]
	for j in range(len(a)):
		a[j] = a[j].split('\t')
		if a[j][3] in MotifTF[0]:
			N[0].append(a[j][1]+'\t'+a[j][3])
			N[1].append(float(a[j][5]))
			N[2].append(MotifTF[1][MotifTF[0].index(a[j][3])])
	VarTF = [[],[]]
	for j in range(len(O[0])):
		if O[0][j] in N[0]:
			indel = N[0].index(O[0][j])
			Fold = np.abs(np.log((O[1][j]+0.01)/(N[1][indel]+0.01)))
			if Fold != 0:
				VarTF[0].append(O[2][j]+'('+str(Fold)+')')
				VarTF[1].append(Fold)
		else:
			Fold = np.abs(np.log((O[1][j]+0.01)/(0.01)))
			VarTF[0].append(O[2][j]+'('+str(Fold)+')')
			VarTF[1].append(Fold)
	for j in range(len(N[0])):
		if N[0][j] not in O[0]:
			Fold = np.abs(np.log((N[1][j]+0.01)/(0.01)))
			VarTF[1].append(Fold)
			VarTF[0].append(N[2][j]+'('+str(Fold)+')')
	if len(VarTF[0]) > 0:
		MaxV[0].append(RE_SNP[i][3])
		MaxV[1].append(max(VarTF[1]))
		SortedVarTF = []
		for j in range(len(VarTF[0])):
			MaxT = max(VarTF[1])
			indel = VarTF[1].index(MaxT)
			SortedVarTF.append(VarTF[0][indel])
			VarTF[1][indel] = -1
		MaxV[2].append(';'.join(SortedVarTF))
if len(MaxV[1]) > 0:
	FinalMaxV = max(MaxV[1]);Final_rs = ''
	for i in range(len(MaxV[0])):
		MaxVV = max(MaxV[1])
		indel1 = MaxV[1].index(MaxVV);
		Final_rs += (str(MaxV[0][indel1])+':'+MaxV[2][indel1]+';;')
		MaxV[1][indel1] = -100
	print(str(FinalMaxV)+'\t'+Final_rs.strip(';;'))
else:
	print('0.0\tNA')
