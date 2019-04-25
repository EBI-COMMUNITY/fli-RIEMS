##### Summerize Results of Blastp 
##### Get variables
arbeitsvz <- as.character(Sys.getenv("arbeitsvz"))
ContigRead <- as.character(Sys.getenv("ContigRead"))

if(ContigRead == "Contig")
{
  workdir <- paste(arbeitsvz,"/assemblyFA/", sep = "")
} else if(ContigRead == "Read") {
  Blastordner <- as.character(Sys.getenv("blastordner"))
  workdir <- Blastordner
}

setwd(workdir)

  Blastp <- read.table("Hits.txt", header = F, stringsAsFactors = F, sep = "\t", dec = ".")
  Blastp$V11 <- Blastp$V1
  Blastp$V1 <- gsub("_([^_]*)$", "", Blastp$V1)
  Blastp$V12 <- Blastp$V9*Blastp$V10/100
  Blastp <- Blastp[ order(Blastp$V1, -Blastp$V12),]
  Blastp <- Blastp[!duplicated(Blastp$V1),]
  Blastp[,12] <- NULL
  Blastp[,11] <- NULL
  Blastp[,10] <- NULL
  Blastp[,9] <- NULL
  
write.table(file = "Hits.txt", x = Blastp, quote = F, row.names = F, col.names = F, sep = "\t")











