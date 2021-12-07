from Bio import SeqIO
import sys
import pandas
import random

fq=sys.argv[1]
output=sys.argv[2]

lst=[]

perc=.01

for record in SeqIO.parse(fq,"fastq"):
	lst.append(record.letter_annotations["phred_quality"])

readLen=max([len(l) for l in lst])

print(readLen)


lst=[l+[-1]*(readLen-len(l)) for l in lst]

cols=["pos"+str(i) for i in range(0,readLen)]

df=pandas.DataFrame(lst, columns =cols)

df.to_csv(output,sep="\t")
