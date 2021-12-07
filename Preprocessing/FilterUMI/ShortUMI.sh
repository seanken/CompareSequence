##Replaces the last 3 bases of the UMI with A's and updates the quality of those based to be Q40
fq1=$1 #Read 1
fq2=$2 # Read 2
outdir=$3 #The output directory

mkdir $outdir

zcat $fq1 | awk '{if(NR%4==2){print substr($0, 1, length($0)-3)"AAA"}else{print$0}}' | awk '{if(NR%4==0){print substr($0, 1, length($0)-3)"III"}else{print$0}}'  | gzip > $outdir/out_S1_L001_R1_001.fastq.gz
ln -s $fq2 $outdir/out_S1_L001_R2_001.fastq.gz
