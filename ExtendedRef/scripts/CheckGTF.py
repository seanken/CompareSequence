import sys

##Cleans GTF to  make it match
def checkGTF(gtf_old):
	print("Read in!")
	fil=open(gtf_old,"r")
	dct={}
	for lin in fil:
		lin=lin.strip()
		s=lin.split("\t")
		if len(s)<8:
			continue;
		tags=s[8]
		chrom=s[0]
		if s[2]!="exon":
			continue;
		gene_id=tags.split("gene_id ")[1].split(" ")[0].strip(";").strip("\"")
		gene_name=tags.split("gene_name ")[1].split(" ")[0].strip(";").strip("\"")
		ID=gene_id+"__"+gene_name
		curChrom=dct.get(ID,chrom)
		if curChrom!=chrom:
			chrom="NOTTA"
		dct[ID]=chrom
	fil.close()
	print(len([i for i in dct]))
	print(len([i for i in dct if dct[i]=="NOTTA"]))


if __name__=="__main__":
	gtf_old=sys.argv[1]
	checkGTF(gtf_old)
