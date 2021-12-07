
#loads all 10X lanes from a given directort
#'Create Seurat object
#'
#'Creates a Seurat object given a count matrix or reprocesses an existing Seurat object
#'
#' @param seur, a preexisting Seurat object
#' @param dat, a count Matrix (genes by cells)
#' @param minGenes, the min number of genes per cell
#' @param regress, things to regress out
#' @param species, the species being looked at 
#' @return A Seurat object
#' @export
dir10X<-function(seur=NULL,dat=NULL,minGenes=500,regress=c("nCount_RNA"),species="human")
{

if(is.null(seur))
{
print("Make object!")
seur<-CreateSeuratObject(dat,"Seurat",min.features=minGenes)#,normalization.method="LogNormalize",scale.factor=1000000)
}

seur<-NormalizeData(seur,normalization.method="LogNormalize",scale.factor=1000000)

print("Get variable genes!")
seur<-FindVariableFeatures(seur)


print("Regress out!")
if(length(regress)>0)
{
seur<-ScaleData(seur,features=seur@assays$RNA@var.features,vars.to.regress=regress)

}
else{
seur<-ScaleData(seur,features=seur@assays$RNA@var.features)
}
print("Run PCA!")
seur<-RunPCA(seur,npcs=60)


print("Return!")

return(seur)


}



