f = open('./SpecificGenes/TissueSpeGene.txt')
SpeTF = f.readlines();f.close()
SpeTF = [SpeTF[i].strip('\n') for i in range(len(SpeTF))]

f = open('GWAS_Tissue_rs_sorted.txt')
a = f.readlines();f.close()
g = open('GWAS_Tissue_rs_sorted_clean.txt','w')
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	if a[i][6] != 'NA' and float(a[i][5]) > 0:
		V = []
		a[i][6] = a[i][6].split(';;')
		for j in range(len(a[i][6])):
			a[i][6][j] = a[i][6][j].split(':')
			a[i][6][j][1] = a[i][6][j][1].split(';')
			TF = [];Fold = []
			for k in range(len(a[i][6][j][1])):
				a[i][6][j][1][k] = a[i][6][j][1][k].split('(')
				fold = float(a[i][6][j][1][k][1].strip(')'))
				a[i][6][j][1][k][0] = a[i][6][j][1][k][0].split(',')
				for l in range(len(a[i][6][j][1][k][0])):
					if a[i][6][j][1][k][0][l] in SpeTF:
						if a[i][6][j][1][k][0][l] not in TF:
							TF.append(a[i][6][j][1][k][0][l])
							Fold.append(fold)
						else:
							indel = TF.index(a[i][6][j][1][k][0][l])
							Fold[indel] = max(fold,Fold[indel])
			if len(TF) > 0:
				tmp = a[i][6][j][0]+':'
				for k in range(len(TF)):
					tmp += (TF[k]+'('+str(Fold[k])+'),')
				V.append(tmp.strip(','))
		if len(V) > 0:
			V = ';'.join(V)
			g.write(a[i][0]+'\t'+a[i][1]+'\t'+a[i][2]+'\t'+a[i][3]+'\t'+a[i][4]+'\t'+a[i][5]+'\t'+V+'\n')
g.close()
