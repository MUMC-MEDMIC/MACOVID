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
   echo "t     number of cores"
   echo
}

######


######
#Options
######

while getopts ":hi::o::m::t::" option; do
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
      t) # threads
         threads=$OPTARG;;
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
python macovid.py mapreads -i $input -o $output -m $mfile --cores $threads --scheme primer_schemes/EMC/v2/