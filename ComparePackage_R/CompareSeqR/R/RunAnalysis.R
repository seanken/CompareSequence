#' Generates Seurat objects for comparision
#'
#' Feed in two sequencing runs on same cells, outputs plots, etc. Note final plots in manuscript might be slightly different (modified later) but the underlying data should be the same.
#'
#' @param x a character array of length 5. The first entry is the sample name, the second the CellRanger directory of seq run 1, the third the CellRanger directory of seq run 2, the 4th the output of vireo on seq run 1 (N if not used), and the 5th the output of vireo on seq run 2 (N if not used)
#' @param outdir The directory to output results to
#' @return lst List of images
#' @export
ProcessCellRangerOuput=function(x,outdir)
{
sampname=x[1]
ill=x[2]
ult=x[3]

vir_ill=x[4]
vir_ult=x[5]


print("Next sample!")
system(paste("mkdir ",outdir,sep=""))

print("Compare CellRanger Output Cells")
bars_ill=scan(paste(ill,"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",sep=""),"")
bars_ult=scan(paste(ult,"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",sep=""),"")

num_ill=length(bars_ill)
num_ult=length(bars_ult)
num_both=length(intersect(bars_ill,bars_ult))

dat=data.frame(c((num_ill-num_both),(num_ult-num_both),num_both))
colnames(dat)[1]="Number_of_Cells"
dat["Type"]=c("Illumina Only","Ultima Only","Both")
p1=ggplot(dat,aes(x=Type,y=Number_of_Cells))+geom_bar(stat="identity")+coord_flip()+ggtitle("Number of shared cells")+xlab("")+ylab("Number of Cells")
write.table(dat,paste(outdir,"/cell.counts.txt",sep=""),row.names=F,quote=F)

print("Get WC info")
dat_ill=Read10X(paste(ill,"/outs/filtered_feature_bc_matrix",sep=""))
dat_ult=Read10X(paste(ult,"/outs/filtered_feature_bc_matrix",sep=""))

inter=intersect(colnames(dat_ill),colnames(dat_ult))
dat_ill=dat_ill[,inter]
dat_ult=dat_ult[,inter]
nGene_ill=colSums(dat_ill>0)
nGene_ult=colSums(dat_ult>0)
nUMI_ill=colSums(dat_ill)
nUMI_ult=colSums(dat_ult)

tab_ill=data.frame(cbind(nGene_ill,nUMI_ill))
tab_ult=data.frame(cbind(nGene_ult,nUMI_ult))
colnames(tab_ill)=c("nGene","nUMI")
colnames(tab_ult)=c("nGene","nUMI")
tab_ill["Type"]="Illumina"
tab_ult["Type"]="Ultima"

tab=rbind(tab_ill,tab_ult)

p2=ggplot(tab,aes(x=Type,y=nGene,fill=Type))+geom_violin(scale="width")+ggtitle("nGene")+theme(legend.position="none")
p3=ggplot(tab,aes(x=Type,y=nUMI,fill=Type))+geom_violin(scale="width")+ggtitle("nUMI")+theme(legend.position="none")

print("Compare per gene")
mn_ill=rowSums(dat_ill)
mn_ult=rowSums(dat_ult)
tab=cbind(mn_ill,mn_ult)
tab=log(tab+1,10)
tab=data.frame(tab)
colnames(tab)=c("Illumina","Ultima")
p7=ggplot(tab,aes(x=Ultima,y=Illumina))+geom_point()+ggtitle("Gene Expression Comparision")+xlab("nUMI in Ultima")+ylab("nUMI in Illumina")

if(nchar(vir_ill)>5 & nchar(vir_ult)>5)
{
print("Use Vireo")
doub_ill=read.table(vir_ill,header=T,stringsAsFactors=F,sep="\t")
doub_ult=read.table(vir_ult,header=T,stringsAsFactors=F,sep="\t")

doub_ill=doub_ill[grep("^donor",doub_ill[,2]),]
doub_ult=doub_ult[grep("^donor",doub_ult[,2]),]

print(head(doub_ill))
print(head(doub_ult))

inter=intersect(doub_ill[,1],doub_ult[,1])

inter=intersect(inter,colnames(dat_ill))
inter=intersect(inter,colnames(dat_ult))

dat_ill=dat_ill[,inter]
dat_ult=dat_ult[,inter]
}

print("Generate Seurat Objects")
seur_ill=dir10X(dat=dat_ill,minGenes=0,regress="nFeature_RNA")
seur_ill=RunUMAP(seur_ill,dims=1:20)
seur_ult=dir10X(dat=dat_ult,minGenes=0,regress="nFeature_RNA")
seur_ult=RunUMAP(seur_ult,dims=1:20)

print("Label Cell Types with Azimuth")
seur_ill=RunAzimuth_pbmc(seur_ill,"celltype.l1")
seur_ult=RunAzimuth_pbmc(seur_ult,"celltype.l1")

p4=UMAPPlot(seur_ill,group.by="predicted.subclass",label=T)+theme(legend.position="none")+ggtitle("Illumina UMAP")
p5=UMAPPlot(seur_ult,group.by="predicted.subclass",label=T)+theme(legend.position="none")+ggtitle("Ultima UMAP")

print("Get CellType Comp Plot")
tab=cbind(seur_ill@meta.data[,"predicted.subclass"],seur_ult@meta.data[,"predicted.subclass"])
tab=data.frame(tab)
colnames(tab)=c("Illumina_Cell_Types","Ultima_Cell_Types")
tab<-tab %>% group_by(Illumina_Cell_Types,Ultima_Cell_Types) %>% summarise(NumCells=length(Ultima_Cell_Types)) %>% as.data.frame()
tab["Match"]=tab[,"Illumina_Cell_Types"]!=tab[,"Ultima_Cell_Types"]
acc=sum(tab[!tab[,"Match"],3])/sum(tab[,3])
p6=ggplot(tab,aes(y=NumCells,axis1=Illumina_Cell_Types,axis2=Ultima_Cell_Types))+geom_alluvium(aes(fill=Match),width = 0)+geom_stratum(width = 1/12, fill = "white")+geom_text(stat = "stratum", aes(label = after_stat(stratum)))+guides(fill = FALSE)+scale_x_continuous(breaks = 1:2, labels = c("Illumina","Ultima"))+theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.line.x=element_blank(),axis.line.y=element_blank())+ylab("")+xlab("")+scale_fill_manual(values=c("grey","red"))+ggtitle(paste("Cell Type composition, ",floor(100*acc),"% Agreement",sep=""))


p=(p1+p7)/(p2+p3)
p=p/p6
p=p/(p4+p5)

lst=list()
lst[[1]]=p1
lst[[2]]=p2
lst[[3]]=p3
lst[[4]]=p4
lst[[5]]=p5
lst[[6]]=p6
lst[[7]]=p7

saveRDS(lst,paste(outdir,"/plots.RDS",sep=""))

ggsave(paste(outdir,"/plots.pdf",sep=""),p,width=14,height=20)
ggsave(paste(outdir,"/",sampname,"plots.pdf",sep=""),p,width=14,height=20)


print("Save Metadata objects!")
saveRDS(seur_ill@meta.data,paste(outdir,"/seur.meta.ill.RDS",sep=""))
saveRDS(seur_ult@meta.data,paste(outdir,"/seur.meta.ult.RDS",sep=""))

print("Save Seurat objects!")
saveRDS(seur_ill,paste(outdir,"/seur.ill.RDS",sep=""))
saveRDS(seur_ult,paste(outdir,"/seur.ult.RDS",sep=""))

return("Yay!")

}
