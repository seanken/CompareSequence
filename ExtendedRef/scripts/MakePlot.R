library(Seurat)
library(ggrepel)
library(Matrix)
library(ggplot2);library(patchwork);library(cowplot)
theme_set(theme_cowplot())


MakePlot=function(dir1,dir2,plot_title="")
{
print("load data")
dat_ill=Read10X(paste(dir1,"/outs/filtered_feature_bc_matrix",sep=""))
dat_ult=Read10X(paste(dir2,"/outs/filtered_feature_bc_matrix",sep=""))
inter=intersect(rownames(dat_ill),rownames(dat_ult))
dat_ill=dat_ill[inter,]
dat_ult=dat_ult[inter,]

inter=intersect(colnames(dat_ill),colnames(dat_ult))

dat_ill=dat_ill[,inter]
dat_ult=dat_ult[,inter]

mn_ill=rowSums(dat_ill)
mn_ult=rowSums(dat_ult)
tab=cbind(mn_ill,mn_ult)
#tab2=tab
#tab=log(tab+1,10)
tab=data.frame(tab)
colnames(tab)=c("Illumina","Ultima")

dat=tab;
dat[1]=1000000*dat[,1]/sum(dat[,1]);dat[2]=1000000*dat[,2]/sum(dat[,2]);dat["Rat"]=abs(log(dat[,1]+10)-log(dat[,2]+10));dat[1]=log(dat[,1]+1);dat[2]=log(dat[,2]+1);dat["Diff"]=dat[,"Rat"]>log(2);
tab=dat
tab["Gene"]=rownames(dat_ill)
print(head(dat))
tab3=tab[tab$Rat>log(2),]
tab3=tab3[order(tab3$Rat,decreasing=T),]
print(head(tab3))

p7=ggplot(tab,aes(x=Ultima,y=Illumina,color=Diff))+geom_point()+ggtitle(plot_title)+xlab("log TPM Standard Reference")+ylab("log TPM Extended Reference")+theme(legend.position="none")+scale_fill_manual(values=c("red","black"))+geom_text_repel(aes(x=Ultima,y=Illumina,label=Gene),data=tab3[1:20,],color="black",max.overlaps=20)

return(p7)

}




