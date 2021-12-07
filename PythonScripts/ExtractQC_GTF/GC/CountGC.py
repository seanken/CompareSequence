fil=open("temp.fa")
dct_gc={}
dct_len={}
for line in fil:
	gene=line.split(":")[0]
	seq=line.split("\t")[1].strip()
	numGC=len([s for s in seq if s in ["G","C"]])
	Len=len([s for s in seq])
	dct_gc[gene]=dct_gc.get(gene,0)+numGC
	dct_len[gene]=dct_len.get(gene,0)+Len
fil.close()

print("save")
savfil=open("GC.counts.txt","w")
i=0
for gene in dct_gc.keys():
	i=i+1
	perc=dct_gc[gene]/dct_len[gene]
	toWrite=gene+"\t"+str(dct_gc[gene])+"\t"+str(dct_len[gene])+"\t"+str(perc)+"\n"
	if i<10:
		print(toWrite)
	savfil.write(toWrite)
savfil.close()
