#!/bin/bash
#SBATCH --job-name=macovid
#SBATCH --output=macovid.%j.txt

######
#Help display
######


Help()
{
   # Display Help
   echo "Add description of the script functions here."
   echo
   echo "options:"
   echo "h     Print this Help."
   echo "i     Input directory basecalled data."
   echo "o     Output directory covid genomes."
   echo "m     manifest file."
   echo "p     scheme prefix."
   echo "s     primer_scheme path."
   echo "t     number of cores."
   echo "v     vcf majority call (0-100)"
   echo
}

######


######
#Options
######

#set defaults
prefix=nCoV-2019
scheme=primer_schemes/EMC/V3
majority=66

while getopts ":hi::o::m::p::s::t::v::" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      i) # input directory
         input=$OPTARG;;
      o) # output directory
         output=$OPTARG;;
      m) # manifest file
         mfile=$OPTARG;;
      p) # prefix
         prefix=$OPTARG;;
      s) # primer_scheme path
         scheme=$OPTARG;;
      t) # threads
         threads=$OPTARG;;
      v) # vcf filter
         majority=$OPTARG;;
      \?) # incorrect option
         echo "Error: Invalid option\n"
         Help
         exit;;
   esac
done
######

source ~/.bashrc

######
#concatenate fastq.gz files
######

echo "concatenate fastq files"
rm ${input}/*.fastq
for i in $(ls -d ${input}/*); do echo $i; zcat ${i}/*.fastq.gz > ${i}.fastq; done
echo "finished"

######
#run macovid
######

conda activate macovid
python macovid.py mapreads -i $input -o $output -m $mfile --cores $threads --scheme $scheme --scheme_prefix $prefix --majority $majority