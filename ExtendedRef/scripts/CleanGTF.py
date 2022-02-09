import sys

##Cleans GTF to  make it match
def cleanGTF(gtf_old,gtf_new):
	print("Read in!")
	fil=open(gtf_old,"r")
	savfil=open(gtf_new,"w")
	dct={}
	for lin in fil:
		lin=lin.strip()
		s=lin.split("\t")
		chrom=s[0]
		tags=s[8]
		if s[2]!="exon":
			continue;
		gene_id=tags.split("gene_id ")[1].split(" ")[0].strip(";").strip("\"")
		gene_name=tags.split("gene_name ")[1].split(" ")[0].strip(";").strip("\"")
		trans_name=tags.split("transcript_name ")[1].split(" ")[0].strip(";").strip("\"")
		trans_id=tags.split("transcript_id ")[1].split(" ")[0].strip(";").strip("\"")
		trans_name=trans_name+"_"+gene_id
		trans_id=trans_id+"_"+gene_id
		ID=gene_id+"__"+gene_name
		dct[ID]=dct.get(ID,0)+1 #the exon number
		num=dct[ID]
		num=str(num)
		s[8]="gene_id \""+gene_id+"\"; gene_name \""+gene_name+"\"; gene_version \"6\"; gene_source \"ensembl_havana\"; gene_biotype \"lincRNA\"; transcript_version \"5\"; transcript_source \"havana\"; transcript_biotype \"lincRNA\"; transcript_id \""+trans_id+"\"; transcript_name \""+trans_name+"\"; exon_id \"Exon_"+gene_id+"_"+num+"\"; exon_number \""+num+"\";"
		lin="\t".join(s)
		lin=lin+"\n"
		savfil.write(lin)
	fil.close()
	savfil.close();


if __name__=="__main__":
	gtf_old=sys.argv[1]
	gtf_new=sys.argv[2]
	cleanGTF(gtf_old,gtf_new)
