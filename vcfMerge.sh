scheme=$1
pool1Vcf=$2
pool2Vcf=$3
curDir=$4
outDir=$5
prefix=$6

cd $outDir
artic_vcf_merge $prefix ${curDir}/${scheme} nCoV-2019_1:${pool1Vcf} nCoV-2019_2:${pool2Vcf}
cd $curDir
bgzip -f ${outDir}${prefix}.merged.vcf
tabix -p vcf ${outDir}${prefix}.merged.vcf.gz