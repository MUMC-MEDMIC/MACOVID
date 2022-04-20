#read all arguments
args <- commandArgs(trailingOnly = TRUE)

medakaVcf <- args[1]
recallVcf <- args[2]
longshotVcf <- args[3]
majority <- as.numeric(args[4])
vcfPass <- args[5]
vcfFail <- args[6]

#read Vcf files
medakaVcf <- readLines(medakaVcf)
recallVcf <- readLines(recallVcf)
longshotVcf <- readLines(longshotVcf)

#extract vcf table
medakaVcf.table <- read.table(text=medakaVcf[(max(grep("#",medakaVcf))+1):length(medakaVcf)],sep = "\t")
recallVcf.table <- read.table(text=recallVcf[(max(grep("#",recallVcf))+1):length(recallVcf)],sep = "\t")
longshotVcf.table <- read.table(text=longshotVcf[(max(grep("#",longshotVcf))+1):length(longshotVcf)],sep = "\t")

#check for large deletions
large.del <- medakaVcf.table[nchar(as.vector(medakaVcf.table[,4]))>50,]

#classify SNPs and INDEL
recallVcf.table[,4] <- as.character(recallVcf.table[,4])
recallVcf.table[,5] <- as.character(recallVcf.table[,5])
recallVcf.table[nchar(recallVcf.table[,4]) == nchar(recallVcf.table[,5]),11] <- "mSNP"
recallVcf.table[nchar(recallVcf.table[,4]) == 1 & nchar(recallVcf.table[,5]) == 1,11] <- "SNP"
recallVcf.table[nchar(recallVcf.table[,4]) > nchar(recallVcf.table[,5]),11] <- "DEL"
recallVcf.table[nchar(recallVcf.table[,4]) < nchar(recallVcf.table[,5]),11] <- "INS"

#concatenate recall and longshot vcf files
if(nrow(recallVcf.table) > 1 & nrow(longshotVcf.table) > 1){
	x <- rbind(longshotVcf.table,recallVcf.table[recallVcf.table[,11] %in% c("INS","DEL"),1:10])
} else {
	x <- recallVcf.table
}

#extract number of reads for majority call
x[,11] <- apply(x, 1, function(x){strsplit(gsub("AC=","",strsplit(as.character(x[8]),split = ";")[[1]][2]),split = ",")[[1]][1]})
x[,12] <- apply(x, 1, function(x){strsplit(gsub("AC=","",strsplit(as.character(x[8]),split = ";")[[1]][2]),split = ",")[[1]][2]})


#calculate ALT/REF reads
x[,13] <- apply(x,1,function(x){as.numeric(x[12]) / (as.numeric(x[11]) + as.numeric(x[12])) * 100})

#variants passing majority call
x[,14] <- as.numeric(x[,13]) >= majority

#split vcf in pass and fail
vcfPass.table <- x[x[,14] %in% TRUE,]
vcfFail.table <- x[x[,14] %in% FALSE,]
vcfPass.table <- vcfPass.table[!(vcfPass.table[,2] %in% vcfFail.table[,2]),]

#remove failed ALT that meet majority call for REF
vcfFail.table <- vcfFail.table[!(as.numeric(vcfFail.table[,13]) <= (100-majority)),]


#change medaka format to longshot format
if(nrow(large.del)>0)
{
  large.del[,6] <- "500.00" 
  large.del[,7] <- "PASS"
  large.del[,8] <- "DP=395;AC=0,371;AM=24;MC=0;MF=0.000;MB=0.000;AQ=14.57;GM=1;PH=6.02,6.02,6.02,6.02;SC=None;"
  large.del[,9] <- "GT:GQ:PS:UG:UQ"
  large.del[,10] <- "1/1:500.00:.:1/1:500.00"
}

#add long gap to vcf file
if(nrow(large.del)>0)
{
  vcfPass.table <- rbind(vcfPass.table[,1:10], large.del)
} 

#order vcf files
vcfPass.table <- vcfPass.table[order(vcfPass.table[,2], -nchar(as.character(vcfPass.table[,4]))),1:10]  
vcfFail.table <- vcfFail.table[order(vcfFail.table[,2], -nchar(as.character(vcfFail.table[,4]))),1:10] 

#Vcf header
longshot.header <- longshotVcf[grep("#",longshotVcf)]

#write vcfPass
writeLines(longshot.header,
           args[5])
if(nrow(vcfPass.table) > 0){
	write.table(vcfPass.table,
				args[5],
				append = TRUE,
				quote = FALSE,
				col.names = FALSE,
				row.names = FALSE,
				sep = "\t")
}

#write vcfFail
writeLines(longshot.header,
           args[6])
if(nrow(vcfFail.table) > 0){
	write.table(vcfFail.table,
				args[6],
				append = TRUE,
				quote = FALSE,
				col.names = FALSE,
				row.names = FALSE,
				sep = "\t")
}