import numpy as np

f = open('./Results/GWAS1.SpecVar.RS.txt')
a = f.readlines();f.close()
context = []
RS1 = []
for j in range(len(a)):
	a[j] = a[j].split('\t')
	a[j][1] = float(a[j][1])
	RS1.append(a[j][1])
	context.append(a[j][0])
Rank1 = [0 for i in range(len(a))]
for j in range(len(a)):
	indel = RS1.index(max(RS1))
	Rank1[indel] = j+1
	RS1[indel] = -100000000000

f = open('./Results/GWAS2.SpecVar.RS.txt')
a = f.readlines();f.close()
RS2 = []
for j in range(len(a)):
	a[j] = a[j].split('\t')
	a[j][1] = float(a[j][1])
	RS2.append(a[j][1])
Rank2 = [0 for i in range(len(a))]
for j in range(len(a)):
	indel = RS2.index(max(RS2))
	Rank2[indel] = j+1
	RS2[indel] = -100000000000

Rank = [(Rank1[i]+Rank2[i])/2 for i in range(len(a))]
print(context[Rank.index(min(Rank))])
