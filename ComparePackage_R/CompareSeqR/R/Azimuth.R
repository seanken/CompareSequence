##
##Runs Azimuth on a seurat object, seur, using the reference Azimuth object (created with GenerateReference from a Seurat object)
##Assumes the reference of interest is in reference_map@meta.data in the subclass column
##
#'Run Azimuth with Reference
#'
#'Applies azimuth to a Seurat object with a given reference reference
#'
#' @param seur the Seurat object to label cell types for 
#' @param reference_map Azimuth reference using the subclass column in the reference objects meta data
#' @return Seurat object with Azimuth results and new assay for SCTransform
#' @export
RunAzimuth<-function(seur,reference_map)
{

print("Apply SCTransform")
seur<-SCTransform(seur,new.assay.name = "refAssay",residual.features = rownames(x = reference_map),reference.SCT.model = reference_map[["refAssay"]]@SCTModel.list$refmodel,method = 'glmGamPoi',ncells = 2000,n_genes = 2000,do.correct.umi = FALSE,do.scale = FALSE,do.center = TRUE)

print("Get anchors")
anchors <- FindTransferAnchors(reference = reference_map,query=seur,k.filter = NA,reference.neighbors = "refdr.annoy.neighbors",reference.assay = "refAssay",query.assay = "refAssay",reference.reduction = "refDR",normalization.method = "SCT",features = intersect(rownames(x = reference_map), VariableFeatures(object = seur)),dims = 1:50,n.trees = 20,mapping.score.k = 100)

print("Move over annotation!")
refdata <- lapply(X = "subclass", function(x) {
  reference_map[[x, drop = TRUE]]
})
names(x = refdata) <- "subclass"
if (FALSE) {
 refdata[["impADT"]] <- GetAssayData(
 object = reference_map[['ADT']],
 slot = 'data')
	      }

seur <- TransferData(reference = reference_map,query=seur,dims = 1:50,anchorset = anchors,refdata = refdata,n.trees = 20,store.weights = TRUE)

return(seur)
}


##
##Runs Azimuth on the given Seurat object, seur, using the built in cortex data
##
#'Run Azimuth on human mouse
#'
#'Applies azimuth to a Seurat object with the built in mouse cortex reference
#'
#' @param seur the Seurat object to label cell types for)
#' @return Seurat object with Azimuth results and new assay for SCTransform
#' @export
RunAzimuth_mm10_cortex<-function(seur)
{
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/mouse_motorcortex")
seur=RunAzimuth(seur,reference$map)
return(seur)
}

##
##Runs Azimuth on the given Seurat object, seur, using the built in human cortex data
#'Run Azimuth on human cortex
#'
#'Applies azimuth to a Seurat object with the built in human cortex reference
#'
#' @param seur the Seurat object to label cell types for)
#' @return Seurat object with Azimuth results and new assay for SCTransform
#' @export
RunAzimuth_human_cortex<-function(seur)
{
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_motorcortex")
seur=RunAzimuth(seur,reference$map)
return(seur)
}

##
##Runs Azimuth on the given Seurat object, seur, using the built in human PBMC data
##
#'Run Azimuth on PBMC
#'
#'Applies azimuth to a Seurat object with the built in human PBMC reference
#'
#' @param seur the Seurat object
#' @param col the column in the meta data to use (celltype.l1, celltype.l2, or celltype.l3 all being difference options)
#' @return Seurat object with Azimuth results and new assay for SCTransform
#' @export
RunAzimuth_pbmc<-function(seur,col="celltype.l2")
{
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")
reference$map@meta.data["subclass"]=reference$map@meta.data[,col]
seur=RunAzimuth(seur,reference$map)
return(seur)
}


##
##generate a azimuth reference from a reference seurat object, ref, where meta.data is the column in ref@meta.data containing the reference of interest
##
#'Generate a Azimuth Reference
#'
#'Turns a standard seurat obeject into a Azimuth reference, uses the column in mata.data as the column of interest
#'
#' @param ref A seurat object 
#' @param meta.data string representing a column in meta data
#' @return A Azimuth reference
#' @export 
GenerateReference<-function(ref,meta.data)
{
ref@meta.data["subclass"]=ref@meta.data[,meta.data]
print("Perform SCT")
ref<-SCTransform(ref,new.assay.name = "refAssay",method = 'glmGamPoi',ncells = 2000,n_genes = 2000,do.correct.umi = FALSE,do.scale = FALSE,do.center = TRUE)
print("Make Reference")
reference_map=AzimuthReference(ref,refDR = "pca",refAssay = "refAssay",plot.metadata ="subclass",metadata = "subclass")
return(reference_map)
}



