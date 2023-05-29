
library(tibble)
library(dplyr)
library(tidyr)

library(dplyr)
library(purrr)
library(ggplot2)

require(Biostrings) 
library(data.table)

writeFastaref<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"G4_refrc"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}


writeFastaalt<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"G4_altrc"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

ncv<-readRDS("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/NCV/NCV_1.RDS")
#ncv1 <- ncv[sample(nrow(ncv),10000),]

ncv1<-ncv
rm(ncv)
ncv1<-ncv1%>% filter(V12==1 | V12==-1 | V16==1 | V16== -1)
ncv1<-setDT(ncv1)
ncv1[, G4_refrc := ifelse((V12 == -1 | V16 == -1) , sapply( G4_ref, function(x) as.character(
reverseComplement(DNAString(x)))), sapply(G4_ref,function(x) as.character(x)))]


ncv1[, G4_altrc:= ifelse(V12 == -1 | V16 == -1, sapply( G4_alt, function(x) as.character(
reverseComplement(DNAString(x)))), sapply(G4_alt,function(x) as.character(x) ))]

#ncv1[, G4_altrc := ifelse(V16 == -1, sapply( G4_alt, function(x) as.character(reverseComplement(
#DNAString(x)))),sapply(G4_alt,function(x) as.character(x) ))]

#ncv1[, G4_refrc := ifelse(V16 == -1, sapply( G4_ref, function(x) as.character(
#reverseComplement(DNAString(x)))), sapply(G4_ref,function(x) as.character(x) ))]

ncv1$name<-gsub("*_G4_id","",ncv1$name)
ncv1<-separate(data = ncv1, col = name, into = c("left", "right"), sep = "_")

ncv1$start<-as.numeric(ncv1$right)-30
ncv1$end<-as.numeric(ncv1$right)+30
ncv1$name<-paste0(ncv1$left,":",ncv1$start,"-",ncv1$end)
ncv2<-ncv1
writeFastaref(ncv2,"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/NCV/NCV_1_ref_g4.fasta")
writeFastaalt(ncv2,"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/NCV/NCV_1_alt_g4.fasta")

print(nrow(ncv2))
print(nrow(ncv1))
rm(ncv1)


