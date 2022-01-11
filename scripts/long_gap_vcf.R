args = commandArgs(trailingOnly = TRUE)

vcf <- readLines(args[1])

longshot <- readLines(args[2])

vcf.table <- read.table(text=vcf[(max(grep("#",vcf))+1):length(vcf)],sep = "\t")
large.del <- vcf.table[nchar(as.vector(vcf.table[,4]))>50,]

if(nrow(large.del)>0)
{
  large.del[,6] <- "500.00" 
  large.del[,7] <- "PASS"
  large.del[,8] <- "DP=395;AC=0,371;AM=24;MC=0;MF=0.000;MB=0.000;AQ=14.57;GM=1;PH=6.02,6.02,6.02,6.02;SC=None;"
  large.del[,9] <- "GT:GQ:PS:UG:UQ"
  large.del[,10] <- "1/1:500.00:.:1/1:500.00"
}

longshot.table <- read.table(text=longshot[(max(grep("#",longshot))+1):length(longshot)],sep = "\t")

if(nrow(large.del)>0)
{
  longshot.table <- rbind(longshot.table, large.del)
  longshot.table <- longshot.table[order(as.numeric(longshot.table[,2])),]
}

longshot.header <- longshot[grep("#",longshot)]

writeLines(longshot.header,
           args[3])
write.table(longshot.table,
            args[3],
            append = TRUE,
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)
