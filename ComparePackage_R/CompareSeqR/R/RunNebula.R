#'Run Nebula
#'
#'Runs Nebula on a Seurat object
#' 
#' @param seur the Seurat object to be tested
#' @param form_fixed The formula for the fixed effects (entries should be columns in the meta.data
#' @param sampleCol The column with the sample of origin information (HTO ID in our paper))
#' @param cpc Parameter used to filter out lowly expressed genes, see Nebula documentation
#' @return The DE results in a format from Nebula
#' @export
RunNebula<-function(seur,form_fixed,sampleCol,cpc=.1)
{
meta=seur@meta.data
dat=seur@assays$RNA@counts
print(dim(dat))
print("Reorder")
dat=dat[,order(meta[,sampleCol])]
meta=meta[order(meta[,sampleCol]),]
print("Run DE!")
df = model.matrix(form_fixed, data=meta)
print(head(df))
print(head(meta[,sampleCol]))

re = nebula(dat,meta[,sampleCol],pred=df,offset=meta$nCount_RNA,cpc=cpc)
return(re)
}
