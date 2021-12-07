import pysam;
import random;

##
##Goes read by read to get coverage info for genes in gtf file, looking only at exons
##
def CoverageExons(gtffile,bamfile,perc=.01):
	print("Get Gene Info!")
	geneInfo=getGeneInfo(gtffile);

	print("Iterate over bam file!")
	samfile = pysam.AlignmentFile(bamfile, "rb")
	dist=[]
	iterat=0
	inside=0
	coding=0
	notgeneInfo=0
	notTag=0
	onlyIntron=0
	used=0
	for seq in samfile:
		iterat=iterat+1
		rand=random.uniform(0,1)
		#if iterat % 10000==0:
		#	print(inside)
		#	print(coding)
		#	print(iterat)
		#	print(used)
		#	print(seq.query_name)
		#	print(" ")
		
		if rand>perc:
			continue;
		used=used+1
		lab="BLANK"
		coding=coding+1
		if not seq.has_tag("GN") or not seq.has_tag("CB") or not seq.has_tag("UB"):
			notTag=notTag+1
			continue;
		if not seq.has_tag("NH"):
			continue;
		num=seq.get_tag("NH")
		if num>1:
			continue;
		gene=seq.get_tag("GN")
		if not gene in geneInfo:
			notgeneInfo=notgeneInfo+1
			continue;
		inside=inside+1
		info=geneInfo[gene]
		#sign=seq.get_tag("SG")
		#sign=seq.get_tag("gs").split(",")[0]
		sign=info.sign
		start=seq.reference_start
		end=seq.reference_end 
		

		startperc=info.posExon(start)
		endperc=info.posExon(end)

		if startperc==endperc:
			onlyIntron=onlyIntron+1
			continue;

		if sign=="+":
			temp=1-startperc
			startperc=1-endperc
			endperc=temp

		
		startIn=info.IntronExon(start)
		if sign=="+":
			startIn=info.IntronExon(end)
		
		dist.append([startperc,endperc,gene,sign,info.lenExon,lab,startIn])
	print(onlyIntron)
	print(notTag)
	print(notgeneInfo)
	return(dist)


##For a given gene contains genomic information, including exonic positions
class GeneInfo:
	def __init__(self,geneName,sign,chrom):
		self.geneName=geneName;
		self.sign=sign;
		self.chrom=chrom
		self.start_exons=[]
		self.end_exons=[]
		self.lenExon=-1
		

	##Tells you if a position is an intron, exon, or outside the gene
	def IntronExon(self,pos):
		if pos<min(self.start_exons):
			return("NOT")
		if pos>max(self.end_exons):
			return("NOT")
		for i in range(0,len(self.start_exons)):
			if self.start_exons[i]<=pos and self.end_exons[i]>=pos:
				return("EXON")
		return("INTRON")



	def addExon(self,start,end):
		self.start_exons.append(start)
		self.end_exons.append(end)

	##Gets the proportion of the way into a gene the position is, counting only exonic bases
	def posExon(self,pos):
		distBefore=0
		for i in range(0,len(self.start_exons)):
			if self.end_exons[i]<pos:
				distBefore=distBefore+self.end_exons[i]-self.start_exons[i]+1
				continue;
			if pos<self.start_exons[i]:
				continue;
			distBefore=distBefore+pos-self.start_exons[i]+1
		perc=1.0*distBefore/(self.lenExon)
		return(perc)
			

	def mergeExons(self):
		startNew=[]
		endNew=[]
		comb=[]
		comb.extend(self.start_exons)
		comb.extend(self.end_exons)
		typ=[1 for x in self.start_exons]
		typ.extend([-1 for x in self.end_exons])
		
		tot=[[comb[i],typ[i]] for i in range(0,len(comb))]

		tot.sort()

		scor=0

		for val in tot:
			pos=val[0]
			side=val[1]
			
			if scor==0:
				startNew.append(pos)

			scor=scor+side

			if scor==0:
				endNew.append(pos)

		self.start_exons=startNew;
		self.end_exons=endNew;
		self.lenExon=sum([endNew[i]-startNew[i]+1 for i in range(0,len(startNew))])
			
	


##
##Parses the GTF file to get gene info, in the form of a dictionary mapping from gene names to GeneInfo objects
##
def getGeneInfo(gtffile):
	fil=open(gtffile);
	
	dct={}

	line=fil.readline()
	iterat=0;
	while len(line)>1:
		s=line.split("\t")
		line=fil.readline()
		iterat=iterat+1
		#if iterat % 10000==0:
		#print(iterat)
		#print(len(dct))
		#print(" ")
		if not len(s)>8:
			continue;
		if not s[2] in ["exon","five_prime_utr","three_prime_utr"]:
			continue;

		start=int(s[3])
		end=int(s[4])
		sign=s[6]
		chrom=s[0]
		geneName=s[8].split(" ")[11].strip(";").strip('"')

		curGene=GeneInfo(geneName,sign,chrom)

		if geneName in dct:
			curGene=dct[geneName]



		curGene.addExon(start,end)
		dct[geneName]=curGene
	print("Merge Exons!")
	for geneName in dct.keys():
		dct[geneName].mergeExons();
	
	return(dct)


import sys

if __name__=="__main__":
	print("Run!")
	args=sys.argv

	bamfile=args[1]
	gtffile=args[2]
	savefile=args[3]	
	out=CoverageExons(gtffile,bamfile)
	print("Save!")
	savFil=open(savefile,"w")
	for val in out:
		res=str(val[0])+" "+str(val[1])+" "+val[2]+" "+val[3]+" "+str(val[4])+" "+val[5]+" "+val[6]+"\n"
		savFil.write(res)
	savFil.close()



