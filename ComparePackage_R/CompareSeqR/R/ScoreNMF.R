##
##Run Tests
##
#' Tests NMF 
#'
#' Fits NMF on one seurat object and tests the accuracy on others
#'
#' @param seur_ill The first Seurat object 
#' @param seur_ult The second Seurat object, on the same cells as the first
#' @param numTrain Number of cells to use for training
#' @param numRep The number of times repeated
#' @param k The number of factors to use
#' @return dat A data frame of score information
#' @export
MakeTest<-function(seur_ill,seur_ult,numTrain,numRep=10,k=15)
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
dat_ill_test=seur_ill@assays$RNA@data[,test]
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
score=TestScore(dat_ill_train,dat_ill_test,seed=seed,k=k)
sc_ill=c(sc_ill,score)
print(score)
score=TestScore(dat_ult_train,dat_ill_test,seed=seed,k=k)
sc_ult=c(sc_ult,score)
print(score)
score=TestScore(dat_pert,dat_ill_test,seed=seed,k=k)
sc_pert=c(sc_pert,score)
print(score)


}

dat=data.frame(Ultima=sc_ult,Illumin=sc_ill,Pert=sc_pert)

return(dat)

}


##
##dat1 is training data, dat2 is test data
TestScore<-function(dat1,dat2,k=15,seed=1)
{
w=RunNMF(dat1,k,seed)
val=GetNMFScore(dat2,w,k)
return(val)
}

##Runs NMF
RunNMF<-function(dat,k=15,seed=1)
{
out=nmf(dat,k)
w=out$w
return(w)
}

##Gets the score for NMF projected on the test data
GetNMFScore<-function(A,w,k=15)
{
print("Project")
h=project(A,w=w)
print("MSE")
val=mse(A=A,w=w,d=rep(1,k),h=h)
return(val)
}
