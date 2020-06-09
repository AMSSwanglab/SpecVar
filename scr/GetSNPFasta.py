import os 
import time

f = open('a_snp.txt')
RE_SNP = f.readlines();f.close()

for i in range(len(RE_SNP)):
	RE_SNP[i] = RE_SNP[i].strip('\n').split('\t')
	RE = RE_SNP[i][8]+'_'+RE_SNP[i][9]+'_'+RE_SNP[i][10]
	f = open(RE+'.fa')
	b = f.readlines()
	g = open(RE+'_v_'+RE_SNP[i][3]+'.fa','w')
	g.write(b[0])
	b[1] = list(b[1])
	for j in range(len(b[1])):
		b[1][j] = b[1][j].upper()
	A=b[1][int(RE_SNP[i][1])-int(RE_SNP[i][9])]
	A1=RE_SNP[i][5]
	A2=RE_SNP[i][6]
	if A == A1:
		b[1][int(RE_SNP[i][1])-int(RE_SNP[i][9])]=A2
	if A == A2:
		if len(A1) == 1:
			b[1][int(RE_SNP[i][1])-int(RE_SNP[i][9])]=A1
		else:
			A1 = list(A1)
			for k in range(len(A1)):
				if A1[k] != A2:
					b[1][int(RE_SNP[i][1])-int(RE_SNP[i][9])]=A1[k]
	if A != A1 and A != A2:
		print(A,A1,A2)
	b[1] = ''.join(b[1])
	g.write(b[1])
	g.close()

