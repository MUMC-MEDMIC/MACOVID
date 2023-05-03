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

#classify SNPs and INDEL medakaVcf
medakaVcf.table[,4] <- as.character(medakaVcf.table[,4])
medakaVcf.table[,5] <- as.character(medakaVcf.table[,5])
medakaVcf.table[nchar(medakaVcf.table[,4]) == nchar(medakaVcf.table[,5]),11] <- "mSNP"
medakaVcf.table[nchar(medakaVcf.table[,4]) == 1 & nchar(medakaVcf.table[,5]) == 1,11] <- "SNP"
medakaVcf.table[nchar(medakaVcf.table[,4]) > nchar(medakaVcf.table[,5]),11] <- "DEL"
medakaVcf.table[nchar(medakaVcf.table[,4]) < nchar(medakaVcf.table[,5]),11] <- "INS"

#classify SNPs and INDEL recall
recallVcf.table[,4] <- as.character(recallVcf.table[,4])
recallVcf.table[,5] <- as.character(recallVcf.table[,5])
recallVcf.table[nchar(recallVcf.table[,4]) == nchar(recallVcf.table[,5]),11] <- "mSNP"
recallVcf.table[nchar(recallVcf.table[,4]) == 1 & nchar(recallVcf.table[,5]) == 1,11] <- "SNP"
recallVcf.table[nchar(recallVcf.table[,4]) > nchar(recallVcf.table[,5]),11] <- "DEL"
recallVcf.table[nchar(recallVcf.table[,4]) < nchar(recallVcf.table[,5]),11] <- "INS"

#classify SNPs and INDEL recall
longshotVcf.table[,4] <- as.character(longshotVcf.table[,4])
longshotVcf.table[,5] <- as.character(longshotVcf.table[,5])
longshotVcf.table[nchar(longshotVcf.table[,4]) == nchar(longshotVcf.table[,5]),11] <- "mSNP"
longshotVcf.table[nchar(longshotVcf.table[,4]) == 1 & nchar(longshotVcf.table[,5]) == 1,11] <- "SNP"
longshotVcf.table[nchar(longshotVcf.table[,4]) > nchar(longshotVcf.table[,5]),11] <- "DEL"
longshotVcf.table[nchar(longshotVcf.table[,4]) < nchar(longshotVcf.table[,5]),11] <- "INS"

#extract number of reads for majority call recall (not possible for INDELS)
#longshot does not consistently call INDELS and Medaka does not provide REF and ALT counts
x <- recallVcf.table[recallVcf.table[,11] %in% c("mSNP","SNP"),1:10]
x[,11] <- apply(x, 1, function(x){strsplit(gsub("AC=","",strsplit(as.character(x[8]),split = ";")[[1]][2]),split = ",")[[1]][1]})
x[,12] <- apply(x, 1, function(x){strsplit(gsub("AC=","",strsplit(as.character(x[8]),split = ";")[[1]][2]),split = ",")[[1]][2]})

#calculate ALT/REF reads
x[,13] <- apply(x,1,function(x){as.numeric(x[12]) / (as.numeric(x[11]) + as.numeric(x[12])) * 100})

#variants passing majority call
x[,14] <- as.numeric(x[,13]) >= majority

#depth candidate position
x[,15] <- apply(x, 1, function(x){strsplit(gsub("DP=","",strsplit(as.character(x[8]),split = ";")[[1]][1]),split = ",")[[1]][1]})
x[,15] <- as.numeric(x[,15])

#remove less than 30x coverage SNPs
x <- x[x[,15] >=30,]

#extract number of reads for majority call recall (not possible for INDELS)
#longshot does not consistently call INDELS and Medaka does not provide REF and ALT counts
x2 <- longshotVcf.table[longshotVcf.table[,11] %in% c("mSNP","SNP"),1:10]
x2[,11] <- apply(x2, 1, function(x2){strsplit(gsub("AC=","",strsplit(as.character(x2[8]),split = ";")[[1]][2]),split = ",")[[1]][1]})
x2[,12] <- apply(x2, 1, function(x2){strsplit(gsub("AC=","",strsplit(as.character(x2[8]),split = ";")[[1]][2]),split = ",")[[1]][2]})

#calculate ALT/REF reads
x2[,13] <- apply(x2,1,function(x2){as.numeric(x2[12]) / (as.numeric(x2[11]) + as.numeric(x2[12])) * 100})

#variants passing majority call
x2[,14] <- as.numeric(x2[,13]) >= majority

#depth candidate position
x2[,15] <- apply(x2, 1, function(x2){strsplit(gsub("DP=","",strsplit(as.character(x2[8]),split = ";")[[1]][1]),split = ",")[[1]][1]})
x2[,15] <- as.numeric(x2[,15])

#remove less than 30x coverage SNPs
x2 <- x2[x2[,15] >=30,]


#merge longshot vcf
x <- rbind(x,x2)


#split vcf in pass and fail
vcfPass.table <- x[x[,14] %in% TRUE,]
vcfFail.table <- x[x[,14] %in% FALSE,]
vcfPass.table <- vcfPass.table[!(vcfPass.table[,2] %in% vcfFail.table[,2]),]

#remove failed ALT that meet majority call for REF
vcfFail.table <- vcfFail.table[!(as.numeric(vcfFail.table[,13]) <= (100-majority)),]

#extract INDELS of medakaPass
indel <- medakaVcf.table[medakaVcf.table[,11] %in% c("INS","DEL"),1:10]

#change medaka format to longshot format
if(nrow(indel)>0)
{
  indel[,6] <- "500" 
  indel[,7] <- "PASS"
  indel[,8] <- "DP=400;AC=0,398;AM=2;MC=0;MF=0.000;MB=0.000;AQ=23.21;GM=1;PH=6.02,6.02,6.02,6.02;SC=None;"
  indel[,9] <- "GT:GQ:DP:PS:UG:UQ"
  indel[,10] <- "1/1:500:400:.:1/1:500.00"
}

#add long gap to vcf file
if(nrow(indel)>0)
{
  vcfPass.table <- rbind(vcfPass.table[,1:10], indel)
} 

#order vcf files
vcfPass.table <- vcfPass.table[order(vcfPass.table[,2], -nchar(as.character(vcfPass.table[,4]))),1:10]  
vcfFail.table <- vcfFail.table[order(vcfFail.table[,2], -nchar(as.character(vcfFail.table[,4]))),1:10] 

#Vcf header
longshot.header <- recallVcf[grep("#",recallVcf)]

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