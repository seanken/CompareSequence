##
##Makes the HTO matrix for demuxEM
##input: raw data directory from 10X
##outfile: out file to save to
##
#'Create HTO Matrix
#'
#'Creates a HTO matrix to be used by DemuxEM using a CellRanger output with HTO data stored in the Custom field
#'
#' @param input The input 10X directory (the raw ouptut directory)
#' @param outfile The file the HTO matrix is to be written to
#' @return Writes a HTO matrix ready for DemuxEM
#' @export
makeHTOMat<-function(input,outfile)
{
print("Read In Data!")
dat=Read10X(input,strip.suffix=T)
print("Reformat")
dat=dat[["Custom"]]
mn=colSums(dat)
print(dim(dat))
dat=dat[,mn>0]
print(dim(dat))
dat=data.frame(as.matrix(dat))
dat["Antibody"]=rownames(dat)
dat=dat[,c(dim(dat)[2],1:(dim(dat)[2]-1))]
print("Save!")
write.table(dat,outfile,sep=",",row.names=F,quote=F)
print("Done!")
}


