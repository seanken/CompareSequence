params.cleanpy="$projectDir/scripts/CleanGTF.py"
params.gtf //the input gtf
params.bamtobed="$projectDir/scripts/MakeBed.SC.py"
params.bedtogtf="$projectDir/scripts/MakeGTF.py"
params.bam //the input bulk bam, must be sorted, single end
params.outdir="output"
params.merge_sh="$projectDir/scripts/merge.sh"

//Annotates bam with gtf/FeatureCount
process AnnotateBam 
{
input:
path bam, stageAs:"input.bam" from params.bam
path gtf, stageAs:"genes.gtf" from params.gtf

output:
path "input.bam.featureCounts.bam" into ann_bam

'''
featureCounts -s 1 -t exon -g gene_id -a genes.gtf -o counts.txt -R BAM input.bam
'''
}


//makes bed from bam
process ToBed
{
input:
path bam, stageAs:"input.bam" from ann_bam
path pycode, stageAs:"MakeBed.py" from params.bamtobed

output:
path "split.bed" into bam_bed

'''
echo Hi2
python MakeBed.py input.bam split.init.bed
sort-bed --max-mem 30G split.init.bed > split.bed
'''

}

//merge bed
process MergeBed
{
input:
path split_bed, stageAs:"split.bed" from bam_bed
path merge_sh, stageAs:"merge.sh" from params.merge_sh

output:
path "merge.bed" into merge_bed

'''
echo Hi
source merge.sh split.bed merge.bed
'''
//bedtools merge -s -c 6,6,7,7 -o distinct,count_distinct,count_distinct,distinct -i split.bed > merge.bed


}

//makes GTF adding exons corresponding to region in merged bed, and all lines in the original gtf
process MakeGTF
{
publishDir "${params.outdir}"

input:
path merge_bed, stageAs:"merge.bed" from merge_bed
path gtf, stageAs:"gene.gtf" from params.gtf
path bedtogtf, stageAs:"MakeGTF.py" from params.bedtogtf

output:
path "updated.gtf" into updated_gtf

'''
python MakeGTF.py merge.bed gene.gtf updated.gtf
'''

}


//Cleans GTF to avoid issues making reference
process CleanGTF
{
publishDir "${params.outdir}"

input:
path "updated.gtf" from updated_gtf
path cleanpy, stageAs:"CleanGTF.py" from params.cleanpy

output:
path "clean.sorted.gtf" into clean_gtf

'''
python CleanGTF.py updated.gtf clean.gtf
bedtools sort -i clean.gtf > clean.sorted.gtf
'''
}


