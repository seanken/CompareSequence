##Script to run cellsnp-lite
tendir=$1 #10X CellRanger output directory
vcf=$2 ##vcf file to use
outdir=$3 ##out directory

bam=$tendir/possorted_genome_bam.bam
barfil=$tendir/filtered_feature_bc_matrix/barcodes.tsv.gz
bars=bars.txt
mkdir $outdir
zcat $barfil > $bars
cellsnp-lite -s $bam -b $bars -O $outdir -R $vcf --minMAF 0.1 --minCOUNT 20
