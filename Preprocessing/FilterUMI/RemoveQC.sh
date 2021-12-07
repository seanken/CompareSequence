##Uses awk to remove quality info from fastq input
fq1=$1 #fastq read 1
fq2=$2 #fastq read 2
outdir=$3 #the output directory


zcat $fq1 |  awk '{if(NR%4==0){print substr($0, 1, length($0)-12)"IIIIIIIIIIII"}else{print$0}}'  | gzip > $outdir/out_S1_L001_R1_001.fastq.gz
ln -s $fq2 $outdir/out_S1_L001_R2_001.fastq.gz
