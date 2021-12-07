##
##Run Tests
##
#' Tests NMF 
#'
#' Fits NMF on one seurat object and tests the accuracy on others using cell loadings
#'
#' @param seur_ill The first Seurat object 
#' @param seur_ult The second Seurat object, on the same cells as the first
#' @param numTrain Number of cells to use for training
#' @param numRep The number of times repeated
#' @param k The number of factors to use
#' @return dat A data frame of score information
#' @export
##
MakeTest_Cells<-function(seur_ill,seur_ult,numTrain,numRep=10,k=15)
{
inter=intersect(names(seur_ult@active.ident),names(seur_ill@active.ident))


sc_ill=c()
sc_ult=c()
sc_pert=c()

for(seed in 1:numRep)
{
print(seed)
print("Divide into training and testing")
train=sample(inter,numTrain)
test=setdiff(inter,train)
dat_ill_test=seur_ill@assays$RNA@data[,train]
dat_ill_train=seur_ill@assays$RNA@data[,train]
dat_ult_train=seur_ult@assays$RNA@data[,train]
mn=rowMeans(dat_ill_train)
dat_ill_train=dat_ill_train[mn>.01,]
dat_ill_test=dat_ill_test[mn>.01,]
dat_ult_train=dat_ult_train[mn>.01,]


print("Permute!")
rows=rownames(dat_ill_train)
cols=colnames(dat_ill_train)
dat_pert=as(apply(as.matrix(dat_ill_train),1,sample), "sparseMatrix")
dat_pert=t(dat_pert)
rownames(dat_pert)=rows
colnames(dat_pert)=cols
print("Test!")
score=TestScore_cell(dat_ill_train,dat_ill_test,seed=seed,k=k)
sc_ill=c(sc_ill,score)
print(score)
score=TestScore_cell(dat_ult_train,dat_ill_test,seed=seed,k=k)
sc_ult=c(sc_ult,score)
print(score)
score=TestScore_cell(dat_pert,dat_ill_test,seed=seed,k=k)
sc_pert=c(sc_pert,score)
print(score)


}

#lst=list()
#lst[["ult"]]=sc_ult
#lst[["ill"]]=sc_ill
#lst[["pert"]]=sc_pert

dat=data.frame(Ultima=sc_ult,Illumin=sc_ill,Pert=sc_pert)

return(dat)

}


##
##dat1 is training data, dat2 is test data
TestScore_cell<-function(dat1,dat2,k=15,seed=1)
{
h=RunNMF_cell(dat1,k,seed)
val=GetNMFScore_cell(dat2,h,k)
return(val)
}

##Runs NMF
RunNMF_cell<-function(dat,k=15,seed=1)
{
out=nmf(dat,k)
h=out$h
return(h)
}

##Gets the score for NMF projected on the test data
GetNMFScore_cell<-function(A,h,k=15)
{
print("Project")
w=project(A,h=h)
print("MSE")
val=mse(A=A,w=w,d=rep(1,k),h=h)
return(val)
}
