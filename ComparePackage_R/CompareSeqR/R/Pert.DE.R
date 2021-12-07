##
##Looks to see if the perturbations in the two Seurat object have the same expression patterns
##
#' Perform DE analysis on 2 perturb-seq datasets
#'
#' Perform DE analysis on 2 perturb-seq datasets on the same set of cells (so different seq runs)
#'
#' @param seur1 The first Seurat obeject, should contain DemuxEM info in the Assign and Demux column
#' @param seur2 The first Seurat obeject, should contain DemuxEM info in the Assign and Demux column 
#' @param pert The columns in the meta.data with perturbation information
#' @param ref_pert The perturbation used as the reference in the DE
#' @return ret A list of DE results for seur1 and seur2
#' @export
TestPert<-function(seur1,seur2,ref_pert="ONE_NON-GENE_SITE_18",pert="feature_call")
{
seur1@meta.data["feature_call"]=seur1@meta.data[,pert]
seur2@meta.data["feature_call"]=seur2@meta.data[,pert]
print("Get Same Cells")
seur1=subset(seur1,Demux=="singlet")
seur2=subset(seur2,Demux=="singlet")
print("Get perts")
seur1=subset(seur1,cells=names(seur1@active.ident)[!is.na(seur1@meta.data$feature_call)])
seur2=subset(seur2,cells=names(seur2@active.ident)[!is.na(seur2@meta.data$feature_call)])
lst=table(seur1@meta.data$feature_call)
perts=names(lst)[lst>10]
print("Number of perts:")
print(length(perts))
seur1=subset(seur1,feature_call %in% perts)


lst=table(seur2@meta.data$feature_call)
perts=names(lst)[lst>10]
seur2=subset(seur2,feature_call %in% perts)
print("Number of cells with perts:")
print(length(seur1@active.ident))
print(length(seur2@active.ident))

seur1@meta.data[,"feature_call"]=factor(seur1@meta.data[,"feature_call"])
seur2@meta.data[,"feature_call"]=factor(seur2@meta.data[,"feature_call"])
seur1@meta.data[,"feature_call"]=relevel(seur1@meta.data[,"feature_call"],ref=ref_pert)
seur2@meta.data[,"feature_call"]=relevel(seur2@meta.data[,"feature_call"],ref=ref_pert)


print("Perform DE")
print(perts)
out1=RunNebula(seur1,~feature_call,"Assign")[[1]]
out2=RunNebula(seur2,~feature_call,"Assign")[[1]]
print("Process")
ret=list()
ret[[1]]=out1
ret[[2]]=out2
out=ret
print("Get pvalue")
clean=lapply(out,function(x){x=x[,c("gene",grep("p_feature_call",colnames(x),value=T))] %>% gather(pert,pval,-gene);x["pert"]=sub("p_feature_call","",x[,"pert"]);x=x[order(x$pval),];x["FDR"]=p.adjust(x$pval,"fdr");return(x)})
ret[[3]]=clean[[1]]
ret[[4]]=clean[[2]]

print("Get SE")
clean=lapply(out,function(x){x=x[,c("gene",grep("se_feature_call",colnames(x),value=T))] %>% gather(pert,se,-gene);x["pert"]=sub("se_feature_call","",x[,"pert"]);return(x)})
ret[[5]]=clean[[1]]
ret[[6]]=clean[[2]]

print("Get LogFC")
clean=lapply(out,function(x){x=x[,c("gene",grep("logFC_feature_call",colnames(x),value=T))] %>% gather(pert,logfc,-gene);x["pert"]=sub("logFC_feature_call","",x[,"pert"]);return(x)})
ret[[7]]=clean[[1]]
ret[[8]]=clean[[2]]

print("Combined")
tab1=inner_join(inner_join(ret[[3]],ret[[5]]),ret[[7]])
tab2=inner_join(inner_join(ret[[4]],ret[[6]]),ret[[8]])
ret[[9]]=tab1
ret[[10]]=tab2
#tab1=tab1[order(tab1$pval),]
#tab2=tab2[order(tab2$pval),]

names(ret)=c("Tot_Ill","Tot_Ult","pval_Ill","pval_Ult","se_Ill","se_Ult","logFC_Ill","logFC_Ult","nice_Ill","nice_Ult")

res=list()
res[["DE_seur1"]]=ret[[9]]
res[["DE_seur2"]]=ret[[10]]

return(res)
}
