#' Compare DE between Seurat objects
#'
#' Given two Seurat objects performs DE between them with Presto and reports the results
#'
#' @param seur1 Seurat object number 1
#' @param seur2 Seurat object number 2
#' @return A dataframe with DE reslts with positive values being upregulated in seur2
#' @export
runDE=function(seur1,seur2)
{
print("Get Data")
dat1=seur1@assays$RNA@data
dat2=seur2@assays$RNA@data
inter=intersect(colnames(dat1),colnames(dat2))

dat1=dat1[,inter]
dat2=dat2[,inter]

colnames(dat1)=sub("^","Ill_",colnames(dat1))
colnames(dat2)=sub("^","Ult_",colnames(dat2))

dat=cbind(dat1,dat2)

mn=rowSums(dat)
dat=dat[mn>100,]

print(dim(dat))

print("Normalize")
dat@x <- dat@x / rep.int(colSums(dat), diff(dat@p))

dat@x <- log(dat@x*10000+1)

n=dim(dat1)[2]
print(class(dat))
y=c(rep("Illumina",n),rep("Ultima",n))

print("Run Presto")
out=wilcoxauc(dat, y)

out=out[order(out$logFC),]
out=out[out$group=="Ultima",]


return(out)

}
