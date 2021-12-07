import random;



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
	#num=0
	while len(line)>1:
		s=line.split("\t")
		#if len(s)<3:
		#	continue;
		#if s[2]!="exon":
		#	continue;
		
		line=fil.readline()
		iterat=iterat+1
		#if iterat % 1000==0:
		#	print(iterat)
		#	print(len(dct))
		#	print(" ")
		
		
		if not len(s)>8:
			continue;
		#if not s[2] in ["exon"]:
		
		if not s[2] in ["exon"]:
			continue;

		start=int(s[3])
		end=int(s[4])
		sign=s[6]
		chrom=s[0]
		geneName=s[8].split(" ")[11].strip(";").strip('"')

		curGene=GeneInfo(geneName,sign,chrom)

		if geneName in dct:
			curGene=dct[geneName]
		else:
			if geneName=="ensembl":
				print(line);
				print(s)


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

	gtffile="../../../ReadLocationsInGenes/Code/genes.gtf"
	savefile="Merged.exons.bed"	
	geneInfo=getGeneInfo(gtffile)
	fil=open(savefile,"w")
	i=0
	for gene in geneInfo.keys():
		i=i+1;
		cur=geneInfo[gene]
		start_exons=cur.start_exons
		end_exons=cur.end_exons
		numEx=len(start_exons)
		chrom=cur.chrom
		sign=cur.sign
		for j in range(0,numEx):
			start=start_exons[j]
			end=end_exons[j]
			toWrite=chrom+"\t"+str(start)+"\t"+str(end)+"\t"+gene+":"+str(j)+"\t.\t"+sign+"\t"+gene+"\n"
			if i<10:
				print(toWrite)
			fil.write(toWrite)
	fil.close()

