##Requires DemuxEM to be installed (can do with conda)
rawh5=$1 ##the raw count h5 file
HTOMatrix=$2 ##the HTO Matrix (can be generated from the R package)
output=$3 ##where to output the results

demuxEM --generate-diagnostic-plots $rawh5 $HTOMatrix $output
