import sys
##Makes dictionary based off gtf
def loadGTF(gtf):
	dct={}
	fil=open(gtf)
	for line in fil:
		if line[0]=="#":
			continue;
		s=line.strip().split("\t")
		if s[2]!="gene":
			continue;
		tags=s[8]
		gene=tags.split(" ")[1].strip(";").strip("\"")
		if gene in dct:
			tags="NOPE"
		dct[gene]=tags
	print(len([i for i in dct if dct[i]=="NOPE"]))
	print(len([i for i in dct if dct[i]!="NOPE"]))
	return(dct)


##Makes bed into gtf
##bed is the bedfile to be made into a gtf
##gtf is the gtf used to annotate the bam
##outfile is the outfile
def MakeNewGTF(bed,gtf,outfile):
	print("Load GTF")
	dct=loadGTF(gtf)
	print("Process bed")
	out=open(outfile,"w")
	bedfil=open(bed,"r")
	iterat=0
	for line in bedfil:
		iterat=iterat+1
		if iterat % 10000==0:
			print(iterat)
		s=line.strip().split("\t");
		if int(s[4])!=1 or int(s[5])!=1:
			continue;
		gene=s[6]
		if not gene in dct:
			continue;
		tags=dct[gene]
		if tags=="NOPE":
			continue;
		tags=tags+" transcript_id \"TransID"+str(iterat)+"\"; transcript_version \"1\"; exon_number \"1\"; transcript_name \"TransName"+str(iterat)+"\"; transcript_source \"havana\"; transcript_biotype \"lincRNA\"; exon_id \"Exon"+str(iterat)+"\"; exon_version \"1\";"
		toSave=s[0]+"\thavana\texon\t"+str(int(s[1])-1)+"\t"+str(int(s[2])+1)+"\t.\t"+s[3]+"\t.\t"+tags+"\n"
		out.write(toSave)
	bedfil.close()

	print("Add Old Line")
	gtffil=open(gtf)
	for lin in gtffil:
		if lin[0]=="#":
			continue;
		out.write(lin)
	out.close()



if __name__=="__main__":
	args=sys.argv
	bed=args[1]
	gtf=args[2]
	outfile=args[3]
	MakeNewGTF(bed,gtf,outfile)
