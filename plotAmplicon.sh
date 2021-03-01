input=$1
prefix=$2
curDir=$3
outDir=$4

cd $outDir
artic_plot_amplicon_depth --primerScheme ${curDir}${input} --sampleID $prefix --outFilePrefix $prefix ${prefix}*.depths
cd $curDir
