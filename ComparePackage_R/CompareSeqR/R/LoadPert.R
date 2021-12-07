##Loads the output of Cellranger into a Seurat object, assumes Gene expression and other data process together
##dir is the out directory from CellRanger
#'Load Perturb-seq
#'
#'Loads CellRanger and DemuxEM results into Seurat for Pertub-seq
#'
#' @param dir The CellRanger directory (out directory)
#' @return seur Seurat object with Perturb-seq data
#' @export
LoadPertData<-function(dir)
{
print("Load data frame")
lst=Read10X(paste(dir,"/filtered_feature_bc_matrix",sep=""))

print("Make Seurat Object")
seur=dir10X(dat=lst[["Gene Expression"]],minGenes=0,regress=c())
seur<-RunUMAP(seur,dims=1:20)
seur<-FindNeighbors(seur,dims=1:20);seur<-FindClusters(seur)

print("Add HTO")
seur[["HTO"]]=CreateAssayObject(counts =lst[["Custom"]][,names(seur@active.ident)])
seur <- NormalizeData(seur, assay = "HTO", normalization.method = "CLR")
seur=HTODemux(seur, assay = "HTO", positive.quantile = 0.99)

print("Add DemuxEM")
demux=read.table(paste(dir,"/DemuxEM.txt",sep=""),sep="\t",stringsAsFactors=F,header=T,row.names=1)
rownames(demux)=demux[,"index"]
seur@meta.data["Assign"]=demux[sub("-1","",names(seur@active.ident)),"assign_nice"]
seur@meta.data["Demux"]=demux[sub("-1","",names(seur@active.ident)),"demux_nice"]

print("Add Crispr")
seur[["CRISPR"]]=CreateAssayObject(counts =lst[["CRISPR Guide Capture"]][,names(seur@active.ident)])
tab=read.table(paste(dir,"/crispr_analysis/protospacer_calls_per_cell.csv",sep=""),header=T,sep=",")
rownames(tab)=tab[,1]
tab=tab[intersect(rownames(tab),names(seur@active.ident)),]
seur@meta.data["num_features"]=0
seur@meta.data[rownames(tab),"num_features"]=tab[,"num_features"]
seur@meta.data["feature_call"]=tab[names(seur@active.ident),"feature_call"]


print("Add ADT")
seur[["ADT"]]=CreateAssayObject(counts =lst[["Antibody Capture"]][,names(seur@active.ident)])
DefaultAssay(seur) <- "ADT"
seur <- NormalizeData(seur, normalization.method = "CLR", margin = 2)
DefaultAssay(seur) <- "RNA"


return(seur)

}
