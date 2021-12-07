##
#'Jointly process 2 Seurat objects
#'
#'Processes 2 Seurat objects (one from Illumina, one from Ultima) together
#'
#' @param seur_ill, a Seurat objects to combine, with predicted.subclass columns in their meta data
#' @param seur_ult, a Seurat objects to combine, with predicted.subclass columns in their meta data
#' @return A joint Seurat object
#' @export
JointProcess=function(seur_ill,seur_ult)
{
print("Prep data")
dat1=seur_ill@assays$RNA@counts
colnames(dat1)=sub("^","ill_",colnames(dat1))
dat2=seur_ult@assays$RNA@counts
colnames(dat2)=sub("^","ult_",colnames(dat2))

inter=intersect(rownames(dat1),rownames(dat2))

dat=cbind(dat1[inter,],dat2[inter,])

meta1=seur_ill@meta.data
rownames(meta1)=colnames(dat1)
meta2=seur_ult@meta.data
rownames(meta2)=colnames(dat2)
meta=rbind(meta1,meta2)

print("Process")
seur=dir10X(dat=dat,minGenes=0,regress="orig.ident")
seur=RunUMAP(seur,dims=1:20)
seur@meta.data["CellType"]=meta[names(seur@active.ident),"predicted.subclass"]

return(seur)

}

