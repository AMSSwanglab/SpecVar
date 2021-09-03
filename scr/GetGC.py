from scipy import stats

f = open('./Results/GWAS1.SpecVar.RS.txt')
a = f.readlines();f.close()
a = [float(a[i].split('\t')[1]) for i in range(len(a))]

f = open('./Results/GWAS2.SpecVar.RS.txt')
b = f.readlines();f.close()
b = [float(b[i].split('\t')[1]) for i in range(len(b))]

results = stats.spearmanr(a,b)

g = open('./Results/GWAS1_GWAS2_Corr.txt','w')
g.write("Genetic correlation: "str(results[0])+'\nP-value: '+str(results[1])+'\n')
g.close()
