f = open('./Input/GWAS.bed')
a = f.readlines();f.close()
N = len(a)

f = open('./Data/RELength.txt')
L = f.readlines();f.close()
for i in range(len(L)):
	L[i] = float(L[i].split('\t')[1])
f = open('GWAS_RE_Overlap.txt')
a = f.readlines();f.close()
FE = [[],[]]
for i in range(len(a)):
	a[i] = a[i].split('\t')
	FE[0].append(a[i][0])
	fe = float(a[i][1])/L[i]/(N/3095677412)
	FE[1].append(fe)
g = open('./Results/GWAS/GWAS_FE.txt','w')
g.write("Contexts\tFE\n")
if max(FE[1]) >= 3:
	while max(FE[1]) >= 3:
		indel = FE[1].index(max(FE[1]))
		g.write(FE[0][indel]+'\t'+str(FE[1][indel])+'\n')
		FE[1][indel] = -100
	g.close()
else:
	indel = FE[1].index(max(FE[1]))
	g.write(FE[0][indel]+'\t'+str(FE[1][indel])+'\n')
	g.close()

f = open('./Results/GWAS/GWAS_FE.txt')
celltype = f.readlines();f.close()
del celltype[0]
for i in range(len(celltype)):
	celltype[i] = celltype[i].strip('\n').split('\t')[0]
	f = open('./Data/Networks/'+celltype[i]+'_network.txt')
	net = f.readlines();f.close();del net[0]
	for j in range(len(net)):
		net[j]= net[j].strip('\n').split('\t')
		net[j][3] = net[j][3].split(';')
	f = open('GWAS_RE_Overlap/GWAS_'+celltype[i]+'.txt')
	g = open('./Results/GWAS/GWAS_'+celltype[i]+'_network.txt','w')
	g.write('SNP\tTF\tTG\tRE\tTRS\n')
	a = f.readlines();f.close()
	for j in range(len(a)):
		a[j] = a[j].strip('\n').split('\t')
		for k in range(len(net)):
			if a[j][7] in net[k][3]:
				g.write(a[j][3]+'\t'+net[k][0]+'\t'+net[k][1]+'\t'+a[j][7]+'\t'+net[k][2]+'\n')
	g.close()
