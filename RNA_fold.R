#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  results_file=as.character(args[1])  
}

library(data.table)

print(file)
#file<-"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/pg4_ref_alleleNONCODING_RNAFOLDinput.fasta"

#file<-"/bio/home/goonerrn/g_quadruplex/annotation_database/redoing_reproduce_results/fullfastafile.fasta"
filename <- sub("\\.[[:alnum:]]+$", "", basename(results_file))


file1= dirname(normalizePath(results_file))
print(file1)
file=paste0(file1,"/",filename,".fasta")

added_name_for_output_file<-"RNAFOLD_allG4COSMIC"
mainDir<-paste0("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/RNAfold/",filename,"/")
SubDir <- paste0(mainDir,"com/")


dir.create(file.path( SubDir ), showWarnings = FALSE,recursive=TRUE)

 setwd(SubDir )
 message("Filtering Sequences...")
#reading fasta files
File0 <- unlist(lapply(file, readLines))
n0 <- length(File0)
#getting second line
ID0 <- File0[seq(1,n0,2)]

Sequence0 <- File0[seq(2,n0,2)]
Sequence0.split <- strsplit(Sequence0, "")
Sequence0.length <- sapply(Sequence0.split, function(x) length(x))
#filtering by length size, if less than 5 or greater than 148, removed
a<-5
b<-148
Index0.filter <- which(Sequence0.length >= a & Sequence0.length <= b)
  ID0.filter <- ID0[Index0.filter]
  Sequence0.filter <- Sequence0[Index0.filter]
  File0.filter.size <- c(rbind(ID0.filter, Sequence0.filter))

#temp_prefix<-"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/pg4_ref_allele_RNAFOLDinput"
temp_prefix<-"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/"
name <- paste0(temp_prefix,filename, "_temp",".fasta")
write(File0.filter.size, file=name)


n0 <- length(File0.filter.size)
  ID0 <- File0.filter.size[seq(1,n0,2)]
  Sequence0 <- File0.filter.size[seq(2,n0,2)]
# command1 <- paste0("RNAfold --MEA -d2 -p --noPS -i ",name)
#  RNAfold1 <- system(command1, intern=TRUE)



#'''
##NEED to insert G quad loop here

########
#'''
# Nloop1 <- strsplit(SecondaryStrc, "\\((?=(\\.+\\)))", perl = TRUE)
# Nloop <- listLen(Nloop1) - 1
# Index0.nloop <- which(Nloop == 1)
# ID0.nloop <- ID0[Index0.nloop]
# Sequence0.nloop <- Sequence0[Index0.nloop]
 #File0.filter.nloop <- c(rbind(ID0.nloop, Sequence0.nloop))
 #name <- paste0(prefix, "_filtered.fa")
 #write(File0.filter.nloop, name)

#'''


SubDir <- paste0(mainDir,"com1_1")


dir.create(file.path(SubDir), showWarnings = FALSE,recursive=TRUE)
setwd(file.path( SubDir))


file<-name

message("Calculating Features...")
    command1 <- paste0("RNAfold --MEA  --noconv -d2 -p --noPS -i ", file)
    RNAfold1 <- system(command1, intern=TRUE)


output2<-paste0(SubDir,"/RNAfold_pred.txt")
output2.1<-paste0(SubDir,"/RNAfold_pred_checked.txt")
command2 <- paste0("RNAfold --noPS  --noconv -i ", file, " > " ,output2)
    system(command2)
    command2.1 <- paste0("perl /bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/Kaviar-160204-Public/rnafold/1_check_query_content.pl ",
output2,"  ",output2.1) 
    system(command2.1)


SubDir1 <- paste0(mainDir,"/com_quad")


dir.create(file.path( SubDir1), showWarnings = FALSE,recursive=TRUE)
setwd(file.path( SubDir1))

 

Command_quad<-paste0("RNAfold --MEA  --noconv --gquad -d2 -p --noPS -i ",file)
  RNAfold_quad <- system(Command_quad, intern=TRUE)

   


output2_quad<-paste0(SubDir1,"/RNAfold_pred.txt")
output2.1_quad<-paste0(SubDir1,"/RNAfold_pred_checked.txt")
command2_quad <- paste0("RNAfold --gquad  --noconv --noPS -i ", file, " > ",output2_quad)
    system(command2_quad)
    command2.1_quad <- paste0("perl /bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/Kaviar-160204-Public/rnafold/1_check_query_content.pl ",
 output2_quad,"   ", output2.1_quad)
    system(command2.1_quad)

n <- length(RNAfold_quad)
 SecondaryStrc_quad1 <- RNAfold_quad[seq(3,n,7)]
 SecondaryStrc_quad <- gsub( " .*$", "", SecondaryStrc_quad1)


#_------------
  print("second working")
    #n <- length(RNAfold1)
    #ID <- RNAfold1[seq(1,n,7)]
    #ID <- gsub('>', '', ID)
    #ID <- ID[StemF.ID.index]
  
    #ID
    n <- length(RNAfold_quad)
    ID <- RNAfold_quad[seq(1,n,7)]
    ID <- gsub('>', '', ID)
    #ID <- ID[StemF.ID.index]

 #Sequence
    Sequence <- RNAfold_quad[seq(2,n,7)]
    #Sequence <- Sequence[StemF.ID.index]

library(stringr)


Length <- nchar(Sequence)
    
    #G+C ratio (GC)
    nG <- str_count(Sequence, "G")
    nC <- str_count(Sequence, "C")
    GC <- (nG + nC) / Length
    
    #G/C ratio (G.Cr)
    G.Cr <- nG / nC
    
    #A+T/G+C ratio (AT.GCr)
    nA <- str_count(Sequence, "A")
    nT <- str_count(Sequence, "T")
    AT.GCr <- (nA + nT) / (nG + nC)
    
    #Base proportions Ratios (Ar, Tr, Gr, Cr)
    Ar <- nA / Length
    Tr <- nT / Length
    Gr <- nG / Length
    Cr <- nC / Length


    #Dinucleotide Ratios (DNr)
    AAr <- str_count(Sequence, "AA")/ Length
    GGr <- str_count(Sequence, "GG")/ Length
    CCr <- str_count(Sequence, "CC")/ Length
    TTr <- str_count(Sequence, "TT")/ Length
    AGr <- str_count(Sequence, "AG")/ Length
    ACr <- str_count(Sequence, "AC")/ Length
    ATr <- str_count(Sequence, "AT")/ Length
    GAr <- str_count(Sequence, "GA")/ Length
    GCr <- str_count(Sequence, "GC")/ Length
    GTr <- str_count(Sequence, "GT")/ Length
    CAr <- str_count(Sequence, "CA")/ Length
    CGr <- str_count(Sequence, "CG")/ Length
    CTr <- str_count(Sequence, "CT")/ Length
    TAr <- str_count(Sequence, "TA")/ Length
    TGr <- str_count(Sequence, "TG")/ Length
    TCr <- str_count(Sequence, "TC")/ Length
    
    

 SecondaryStrc1 <- RNAfold_quad[seq(3,n,7)]
    #SecondaryStrc1 <- SecondaryStrc1[StemF.ID.index]
    SecondaryStrc <- gsub( " .*$", "", SecondaryStrc1)
    
    BP <- str_count(SecondaryStrc, "\\(")

    #Pairing Probabilities Structure (PairProbStrc)
    PairProbStrc1 <- RNAfold_quad[seq(4,n,7)]
    #PairProbStrc1 <- PairProbStrc1[StemF.ID.index]
    PairProbStrc <- gsub( " .*$", "", PairProbStrc1)
    
    #RNAfold centroid structure (CentroidStrc)
    CentroidStrc1 <- RNAfold_quad[seq(5,n,7)]
    #CentroidStrc1 <- CentroidStrc1[StemF.ID.index]
    CentroidStrc <- gsub( " .*$", "", CentroidStrc1)
    
    #Maximum Expected Accuracy Structure (MEAStrc)
    MEAStrc1 <- RNAfold_quad[seq(6,n,7)]
   # MEAStrc1 <- MEAStrc1[StemF.ID.index]
    MEAStrc <- gsub( " .*$", "", MEAStrc1)




##################################################
    ######### Secondary Structure Statistics #########
    ##################################################
    
    
    
    #Minimum Free Energy (MFE)
    MFE <- as.numeric(regmatches(SecondaryStrc1,
                                 gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",
                                          SecondaryStrc1, perl=TRUE)))
    
    #Ensemble Free Energy (EFE)
    EFE <- as.numeric(regmatches(PairProbStrc1,
                                 gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",
                                          PairProbStrc1, perl=TRUE)))
    
    #Centroid Free Energy (CFE)
    CFE1 <- sub(".*? (.+)", "\\1", CentroidStrc1)
    CFE1 <- gsub("\\{", "", CFE1)
    CFE1 <- gsub("\\}", "", CFE1)
    CFE1 <- gsub(" ", "", CFE1)
    CFE1 <- unlist(strsplit(CFE1, "d="))
    n2 <- length(CFE1)
    CFE <- as.numeric(CFE1[seq(1,n2,2)])
    
    #Centroid Distance to Ensemble (CDE)
    CDE <- as.numeric(CFE1[seq(2,n2,2)])
    


#Maximum Expected Accuracy Structure Free Energy (MEAFE)

MEAFE1 <- sub(".*? (.+)", "\\1", MEAStrc1)
    MEAFE1 <- gsub("\\{", "", MEAFE1)
    MEAFE1 <- gsub("\\}", "", MEAFE1)
    MEAFE1 <- gsub(" ", "", MEAFE1)
    MEAFE1 <- unlist(strsplit(MEAFE1, "MEA="))
    MEAFE <- as.numeric(MEAFE1[seq(1,n2,2)])
    
    #Maximum Expected Accurary (MEA)
    MEA <- as.numeric(MEAFE1[seq(2,n2,2)])
    
    #Base pair Propensity (BPP)
    BPP <- (BP / Length)
    
    #Frequency of the MFE in ensemble (EFreq)
    EFreq1 <- RNAfold_quad[seq(7,n,7)]
    #EFreq1 <- EFreq1[StemF.ID.index]
    EFreq1 <- unlist(strsplit(EFreq1, ";"))
    EFreq2 <- EFreq1[seq(1,n2,2)]
    
ED1 <- EFreq1[seq(2,n2,2)]
    ED <- as.numeric(regmatches(ED1,
                                gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*", 
                                         ED1, perl=TRUE)))



#Ensemble Diversity (ED)
    ED1 <- EFreq1[seq(2,n2,2)]
    ED <- as.numeric(regmatches(ED1,
                                gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*", 
                                         ED1, perl=TRUE)))
    #Adjusted MFE (MFEadj)
    MFEadj <- (MFE / Length)
    
    #Adjusted EFE (EFEadj)
    EFEadj <- (EFE / Length)
    
    #Adjusted base pair Distance (Dadj)
    Dadj <- (CDE / Length)

SE <- -((Ar*log2(Ar))+(Tr*log2(Tr))+(Gr*log2(Gr))+(Cr*log2(Cr)))
    SEadj <- (SE /Length)
    
    #Difference between MFE and EFE Adjusted (DiffMFE.EFE)
    DiffMFE.EFE <- ((MFE - EFE) / Length)
    
    #Ratio between Adjusted MFE and GC (MFEadj.GC)
    MFEadj.GC <- (MFEadj / GC)
    
    #Ratio between Adjusted MFE and base pairs (MFEadj.BP)
    MFEadj.BP <- (MFEadj / BP)
    
    #Adjusted MEAFE
    MEAFEadj <- (MEAFE / Length)
    
    #Adjusted ED
    EDadj <- (ED / Length)




table <- as.data.frame(cbind(ID,GC,BPP, Dadj, SEadj,
                                 DiffMFE.EFE, Length, AT.GCr,
                                 AAr, GGr, CCr, TTr, AGr, ACr, ATr,
				GAr, GCr, GTr, CAr, CGr, CTr, TAr, TGr, TCr,
                                 MFE,
                                 MFEadj, EFE, EFEadj, CFE, CDE, MEAFE, MEA,
                                 ED, MFEadj.GC, MFEadj.BP, MEAFEadj, EDadj))
    

   #rownames(table) <- ID
   # table[is.na(table)] <- 0
table$name<-ID
 prefix<-"fasta_var"
    final_table <- paste0(mainDir,"/",filename,added_name_for_output_file, prefix,"_features.txt")
    fwrite(table, final_table, sep="\t", quote=F)
    













quad.Structural.Pscore <- function(file, prefix, iterate=1000, threshold=0.1, filter=TRUE){


  MFE <- as.numeric(regmatches(SecondaryStrc1,
                               gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",
                                        SecondaryStrc1, perl=TRUE)))
  iterate=1000
  #MFE P-value (MFE.Pval)
  command3 <- paste0("fasta_ushuffle -n ", iterate, " -k 2 -s 1234 < ", 
                     file, " >", prefix, "_iterated.fa")
  seqs <- system(command3, intern=TRUE)
  
  command4 <- paste0("RNAfold --MEA -d2 -p --noPS -i ", prefix, "_iterated.fa")
  RNAfold2 <- system(command4, intern=TRUE)
  n4 <- length(RNAfold2)
  MFEit1 <- RNAfold2[seq(3,n4,7)]
  MFEit <- as.numeric(regmatches(MFEit1,
                                 gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*", 
                                          MFEit1, perl=TRUE)))
  N <- iterate
  MFEit.along <- seq_along(MFEit)
  MFEit.split <- split(MFEit, ceiling(MFEit.along/N))
  R <- as.vector(sapply(Map(`<=`, MFEit.split, MFE), sum))
  MFE.Pval <- (R / (N + 1))


  MFE.Pval <- round(MFE.Pval, 6)
  table <- as.data.frame(MFE.Pval)
  rownames(table) <- ID
  
  final_table <- paste0(outdir, prefix, "_Structural_Pscore.txt")
  write.table(table, final_table, quote=F, row.names=T, col.names=NA, sep="\t")


#if you want to filter table then
 if(filter == TRUE){
    
    index_thres0 <- which(table$MFE.Pval < threshhold)
    index_thres1 <- (index_thres0*2)-1
    index_thres2 <- index_thres1+1
    index_thres <- c(matrix(c(index_thres1, index_thres2), 2, byrow = T)) 
    File0 <- unlist(lapply(file, readLines))
    File_filter <- File0[index_thres]
    name <- paste0(prefix, "_filtered_Pscore.fa")
    write(File_filter, name)
    
    return(table)
    
  } else if(filter == FALSE) {
    
    return(table)
  }
}
  

