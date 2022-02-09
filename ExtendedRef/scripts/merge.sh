bed_in=$1
bed_out=$2

echo Cluster 
bedtools cluster -s -i $bed_in > clust.bed
awk '{print $1"\t1\t10\t"$7"_"$8"\t"$5"\t"$6"\t"$2"\t"$3"\t"$7}' clust.bed > clust.temp.bed
echo Sort By Cluster
sort-bed --max-mem 30G clust.temp.bed > clust.sort.bed
echo Group by
bedtools groupby -i clust.sort.bed -g 1,4,6,9 -c 8 -o max > end.groupby.bed
bedtools groupby -i clust.sort.bed -g 1,4,6,9 -c 7 -o min > start.groupby.bed
paste start.groupby.bed end.groupby.bed | awk '{print $1"\t"$5"\t"$10"\t"$2"\t.\t"$3"\t"$4}' > comb.bed
echo Sort and reorder
sort-bed --max-mem 30G comb.bed > comb.sort.bed
awk '{print $1"\t"$2"\t"$3"\t"$6"\t1\t1\t"$7}' comb.sort.bed > $bed_out
echo Clean up
#rm clust.bed
#rm clust.temp.bed
#rm clust.sort.bed
#rm end.groupby.bed
#rm start.groupby.bed
#rm comb.bed
#rm comb.sort.bed
