CELL_DATA=$1 #output from cellsnp-lite
outdir=$2 #output directory
n_donor=8

vireo -c $CELL_DATA -N $n_donor -o $outdir
