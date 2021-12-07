import pandas as pd
import sklearn
from sklearn import linear_model



def runMimosca(matfile,metafile,outfile):
	print("Load data")
	mat=pd.read_csv(matfile)
	meta=pd.read_csv(metafile)
	
	print("Run Selection")
	lm=sklearn.linear_model.ElasticNet(l1_ratio=0.5,alpha=0.0005,max_iter=10000)
	#lm=sklearn.linear_model.Lasso()	
	lm.fit(meta,mat)
	B=pd.DataFrame(lm.coef_)
	B.columns=meta.columns
	B.index=mat.columns
	
	

	print("Save File")
	B.to_csv(outfile)


if __name__=="__main__":
	print("Ultima!")
	runMimosca("for.Mimosca.ult.mat.txt","for.Mimosca.ult.meta.txt","B.mim.ult.txt")
	print("Illumina!")
	runMimosca("for.Mimosca.ill.mat.txt","for.Mimosca.ill.meta.txt","B.mim.ill.txt")
	print("Background!")
	for i in range(1,101):
		print(i)
		I=str(i)
		print("Ultima!")
		runMimosca("iterat_ult/iter."+I+".txt.mat.txt","iterat_ult/iter."+I+".txt.meta.txt","iterat_ult/B.EN.ult."+I+".txt")
		print("Illumina!")	
		runMimosca("iterat_ill/iter."+I+".txt.mat.txt","iterat_ill/iter."+I+".txt.meta.txt","iterat_ill/B.EN.ill."+I+".txt")
