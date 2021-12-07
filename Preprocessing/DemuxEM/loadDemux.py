##Code to take output from CellBender and make it useable for R
##only one input, the path to the zarr file
##one output, named as the zarrfile with .txt added
import pandas
import zarr
import sys


zarrFile=sys.argv[1]
outfile=zarrFile+".txt"


print(zarrFile)
print("Load Zarr file")
z2 = zarr.open(zarrFile, mode='r')


print("Make pandas")
df=pandas.DataFrame({'index':z2['GRCh38-rna/barcode_metadata/_index'],'assign':z2['GRCh38-rna/barcode_metadata/assignment'],'type':z2['GRCh38-rna/barcode_metadata/demux_type']})

print("Remove those with demux_type 0")
df=df[df['type']!=0]

print("Rename assingment")
lst=z2['GRCh38-rna/barcode_metadata/_categories/assignment']
df["assign_nice"]=[lst[i] for i in df['assign']]

print("Rename Type")
lst=z2['GRCh38-rna/barcode_metadata/_categories/demux_type']
df["demux_nice"]=[lst[i] for i in df['type']]

print("Save")
df.to_csv(outfile,sep="\t")

print("Done")
