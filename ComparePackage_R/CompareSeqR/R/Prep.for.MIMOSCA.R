library(Seurat)
library(Matrix)

#'Prep for MIMOSCA
#'
#' For a given Seurat object prepares input for the MIMOSCA tool
#'
#' @param seur The input Seurat object. Needs Demux and Assign information from DemuxEM, feature_call to be perturbations in the mata dat
#' @param outfile The output file
#' @return lst a listwhose first item is the matrix, the second the meta data
#' @export
prepData<-function(seur,outfile)
{
print("Cluster")
seur=FindClusters(seur,resolution=.2)
print("remov doublets")
seur=subset(seur,Demux=="singlet")
print("Get perts")
tab=table(seur@meta.data$feature_call)
perts=names(tab)[tab>10]
print(perts)
cells=names(seur@active.ident)[seur@meta.data$feature_call %in% perts]

seur=subset(seur,cells=cells)
print("Make Matrix!")
mat=seur@assays$RNA@data
mn=rowMeans(mat>0)
mat=mat[mn>.05,]
mat=t(mat)
mat=as.matrix(mat)

print("Get Meta!")
meta=seur@meta.data
meta["feature_call"]=factor(meta[,"feature_call"])
#meta["feature_call"]=relevel(meta[,"feature_call"],ref="NO_SITE_1")
meta["nFeature_RNA"]=scale(meta[,"nFeature_RNA"])
meta=model.matrix(~0+nFeature_RNA+feature_call+seurat_clusters+Assign,meta)

print("Save")
write.table(mat,paste(outfile,".mat.txt",sep=""),sep=",",quote=F)
write.table(meta,paste(outfile,".meta.txt",sep=""),sep=",",quote=F)

print("Return")
lst=list()
lst[["Mat"]]=mat
lst[["Meta"]]=meta
return(lst)

}


#'Preps background to process with MIMOSCA
#'
#' preps background input for MIMOSCA
#'
#' @param seur The input Seurat object, see prepData for details
#' @param outdir The directory to output to
#' @param numIter the number of iterations to use
#' @return NA
#' @export
iteratePrep<-function(seur,outdir,numIter=100)
{
system(paste("mkdir",outdir))
for(i in 1:numIter)
{
print(i)
donors=unique(seur@meta.data[,"Assign"])
for(j in donors)
{
seur@meta.data[seur@meta.data[,"Assign"]==j,"feature_call"]=sample(seur@meta.data[seur@meta.data[,"Assign"]==j,"feature_call"])
}
outfil=paste(outdir,"/iter.",i,".txt",sep="")
prepData(seur,outfil)
}
print("Done!")
}

