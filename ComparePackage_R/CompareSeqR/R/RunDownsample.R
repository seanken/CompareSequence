##
##For molinfo, downsamples to different number of reads
##
#'Downsample reads
#'
#' Code to downsample reads to various levels and get number of UMI present at each level
#'
#' @param molinfo The molecular info from CellRanger
#' @param pvals The proportions to try downsampling to
#' @return A list with one entry per proportion with the number of UMIs you get at that level of downsampling
#' @export
dowsampleCurve<-function(molinfo,pvals=.05*1:20)
{
print("Downsample!")
out=lapply(pvals,function(x){
print("Current value:")
print(x)
dat=downsampleReads(molinfo, prop=x)
nUMI=sum(dat)
return(nUMI)
})

print("Clean up!")
out=as.numeric(out)
names(out)=pvals
return(out)
}

