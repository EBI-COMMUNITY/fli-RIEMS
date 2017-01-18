############################### Used functions
options(digits = 3)
arbeitsvz <- as.character(Sys.getenv("arbeitsvz"))
setwd(arbeitsvz)
# is.int()
is.int <- function(x){
  is.integer(x)&&length(x)==0L #gibt TRUE aus wenn x=integer(0)
}

############################### Datasets to analyse

TrimStatus <- read.delim("TrimStatus.txt", quote="", stringsAsFactors=FALSE)
asignedAcc <- read.delim("ViRIEMS_asigned.txt", header=FALSE, quote="", stringsAsFactors=FALSE,col.names = c("Method","SkTax","FamTax","Tax","Scientific.Name","Identity","Accession","ReadStatus","MapIdentity","MapPercent","MapReference"))
FamTaxNames <- read.delim("ViRIEMS_famtax_names.txt",header=FALSE, quote="", stringsAsFactors=FALSE)


############################### Complete Resultprotocol

TrimStatus$SeqId <- TrimStatus$Orig.Trimpoints <- TrimStatus$Orig.Trimmed.Length <- TrimStatus$Used.Trimmed.Length <- NULL
TrimStatus$Trim.Start <- gsub("-.*", "", TrimStatus$Trimpoints.Used)
TrimStatus$Trim.Stop <- gsub(".*-", "", TrimStatus$Trimpoints.Used)
TrimStatus$Trimpoints.Used <- NULL

resultprotocol <- merge(x = TrimStatus,y = asignedAcc,by.x = "Accno",by.y = "Accession",all.x = FALSE,all.y = TRUE)
colnames(resultprotocol) <- c("Accession","Raw.Length","Trimpoint.Start","Trimpoint.End","Method","SkTax","FamTax","Tax","Species","Identity","ReadStatus","MapIdentity","MapPercent","MapReference")
resultprotocol$Identity <- as.numeric(round(resultprotocol$Identity, digits = 2))
#resultprotocol$Identity <- gsub(" ", "",resultprotocol$Identity)

write.table(x = resultprotocol,file = paste(arbeitsvz,"/ViRIEMS_report.txt",sep=""),quote = F,sep = "\t",dec = ".",row.names = F,na = "")

############################### Histogram Protocol



############################### Compact Resultprotocol

resultprotocol$Method <- gsub("-","_",resultprotocol$Method)

TaxIDs <- unique(resultprotocol[order(as.numeric(resultprotocol$SkTax),as.numeric(resultprotocol$FamTax),as.numeric(resultprotocol$Tax)),6:9])
TaxIDs$counts <- 0
unique.methods <- unique(resultprotocol$Method)
TaxIDs$Ident.ViRIEMS <- TaxIDs$ViRIEMS <- 0

#for(i in 1:nrow(TaxIDs)){
#	
#	count <- length(which(TaxIDs$Tax[i]==resultprotocol$Tax))
#	TaxIDs$counts[i] <- count
#}


for(i in 1:nrow(TaxIDs)){
  count <- length(which(TaxIDs$Tax[i]==resultprotocol$Tax))
  TaxIDs$counts[i] <- count
  
  for(j in 1:length(unique.methods)){
    count1 <- length(which(TaxIDs$Tax[i]==resultprotocol$Tax & resultprotocol$Method==unique.methods[j]))
    if(nchar(unique.methods[j])!=0){
      names <- which(colnames(TaxIDs)==unique.methods[j])
      TaxIDs[i,names] <- count1
    }
    blast <- grep("ViRIEMS",unique.methods[j])
    if(is.int(blast)==FALSE & count1!=0){
      
      min.ident <- min(resultprotocol$Identity[which(TaxIDs$Tax[i]==resultprotocol$Tax & resultprotocol$Method==unique.methods[j])],na.rm = T)
      max.ident <- max(resultprotocol$Identity[which(TaxIDs$Tax[i]==resultprotocol$Tax & resultprotocol$Method==unique.methods[j])],na.rm = T)
      
      names2 <- paste("Ident.",unique.methods[j],sep="")
      names3 <- which(colnames(TaxIDs)==names2)
      if(min.ident==max.ident){
        TaxIDs[i,names3] <- paste(min.ident)
      }else{
        TaxIDs[i,names3] <- paste(min.ident,"-",max.ident,sep="")
      }
      
    }
  }
}

TaxIDs <- TaxIDs[,c(4,1:3,5:7)]
TaxIDs <- TaxIDs[!is.na(TaxIDs$SkTax),]
TaxIDs <- TaxIDs[order(TaxIDs$SkTax,TaxIDs$FamTax,TaxIDs$counts,decreasing = T),]

sk <- unique(TaxIDs$SkTax)

Famsnames <- NULL
for(i in 1:length(FamTaxNames[,2])){
  
  if(is.na(FamTaxNames[i,1])==TRUE){
    
    for(i in sk){
      
      j <- which(is.na(TaxIDs$FamTax)==TRUE & TaxIDs$SkTax==i)
      if(is.int(j)==FALSE){
        Famsnames[j] <- c("Unknown",rep(NA,(length(j)-1)))
      }
      
    }
    
  }else{
    j <- which(TaxIDs$FamTax == FamTaxNames[i,1])
    
    Famsnames[j] <- c(FamTaxNames[i,2],rep(NA,(length(j)-1)))
  }
  
  
}

TaxIDs <- cbind(Famsnames,TaxIDs)

colnames(TaxIDs)[1] <- "Family"

write.table(x = TaxIDs,file = "ViRIEMS_report-compact.txt",quote = F,sep = "\t",dec = ".",row.names = F,na = "")

