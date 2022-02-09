import pysam
import sys

def ProcessBam(bamfile,outfile):
	print("Process Bam!")
	samfile = pysam.AlignmentFile(bamfile, "rb")
	iterat=0
	fil=open(outfile,"w")
	for seq in samfile:
		iterat=iterat+1;
		if iterat % 100000==0:
			print(iterat)
		if not seq.has_tag("NH") or not seq.has_tag("XS") or not seq.has_tag("XT"):
			continue;
		num=seq.get_tag("NH")
		if num!=1:
			continue;
		is_assign=seq.get_tag("XS")
		if is_assign!="Assigned":
			continue;
		gene=seq.get_tag("XT")
		
		start=seq.reference_start
		end=seq.reference_end 
		chrom=seq.reference_name
		strand="+" ##strand of gene, assumes R1 is antisense, R2 is sense
		if seq.is_reverse:
			strand="-"
		cigar=seq.cigartuples
		if sum([c[1] for c in cigar if c[0]==3])>10:
			continue;

		toSave=chrom+"\t"+str(start)+"\t"+str(end)+"\tEntry_"+str(iterat)+"\t.\t"+strand+"\t"+gene+"\n"
		fil.write(toSave)
	fil.close()
	samfile.close()


if __name__=="__main__":
	args=sys.argv
	bamfile=args[1]
	outfile=args[2]
	ProcessBam(bamfile,outfile)
		
		

	
