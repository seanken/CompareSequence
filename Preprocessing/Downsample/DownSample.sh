##Downsample reads, assumes reads are named out_S1_L001_R1_001.fastq.gz and out_S1_L001_R2_001.fastq.gz
indir=$1 #indirectory for fastq
outdir=$2 #outdirecotory for fastq
rat=$3 #The ratio (between 0 and 1)


echo $indir
echo $outdir
echo $rat

mkdir $outdir
echo Downsample

fastq1=$indir/out_S1_L001_R1_001.fastq.gz
fastq2=$indir/out_S1_L001_R2_001.fastq.gz
fatq1_out=$outdir/out_S1_L001_R1_001.fastq.gz
fatq2_out=$outdir/out_S1_L001_R2_001.fastq.gz

seqtk sample -s 100 $fastq1 $rat | gzip > $fatq1_out
seqtk sample -s 100 $fastq2 $rat | gzip > $fatq2_out
