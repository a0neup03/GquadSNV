
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



##this is recode old file where the ref and alt are. so we merge them here
codingseq<-readRDS("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/coding/coding_indelfiltered.recode.rds")
NCVseq<-readRDS("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/NCV/NCV_indelfiltered.recode.rds")

NCV<-readRDS("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/NCV/NCV_1.RDS")
coding<-readRDS("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/coding/coding_1.RDS")

###result1
count_coding<-coding[, .(G4hunter_ref = mean(V9), G4hunter_alt = mean(V13),count=.N,sd=sd(V9),sd2=sd(V13)), by = .(V12,V16)]
count_NCV<-NCV[, .(G4hunter_ref = mean(V9), G4hunter_alt = mean(V13),count=.N,sd=sd(V9),sd2=sd(V13)), by = .(V12,V16)]

fwrite(count_coding,sep="\t",file="count_coding.tsv")
fwrite(count_NCV,sep="\t",file="count_NCV.tsv")

temp<-rbind(NCV,coding)
count_temp<-temp[, .(G4hunter_ref = mean(V9), G4hunter_alt = mean(V13),count=.N,sd=sd(V9),sd2=sd(V13)), by = .(V12,V16)]
fwrite(count_temp,sep="\t",file="count_totalSNV.tsv")
rm(temp)
coding_total<-left_join(coding,codingseq,by=c("name","G4_ref","G4_alt","G4_id"))

coding1<-coding_total
rm(coding)
coding1<-coding1%>% filter(V12==1 | V12==-1 | V16==1 | V16== -1)
coding1<-setDT(coding1)
coding1[, G4_refrc := ifelse((V12 == -1 | V16 == -1) , sapply( G4_ref, function(x) as.character(
reverseComplement(DNAString(x)))), sapply(G4_ref,function(x) as.character(x)))]


coding1[, G4_altrc:= ifelse(V12 == -1 | V16 == -1, sapply( G4_alt, function(x) as.character(
reverseComplement(DNAString(x)))), sapply(G4_alt,function(x) as.character(x) ))]
coding1$nameimp<-coding1$name
coding1$name<-gsub("_G4_id","",coding1$name)
coding1<-separate(data = coding1, col = name, into = c("left", "right"), sep = "_")

#did not work
#bedops --ec --element-of 1 <(vcf2bed < ./af-only-gnomad.hg38.vcf.gz ) <(sort-bed  /bio/home/goonerrn/g_quadruplex/annotation_database/KVIAR/COSMIC/seqVariants_cosmic/cosmic_G4_variants.bed) > answer.bed


coding1$start<-as.numeric(coding1$right)-30
coding1$end<-as.numeric(coding1$right)+30

coding2<-coding1
coding2$name<-paste0(coding2$left,":",coding2$start,"-",coding2$end,"||","loc:",coding2$loc_start,"||",coding2$REF,"&",coding2$ALT)


#get variants across two strands
table(coding2$V12,coding2$V16)

rm(coding1,coding,codingseq)

##here we repeat what we did for coding and codingseq with NCV and NCVseq

ref_RNAFOLD_ncv<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/RNAfold/NCV_1_ref_g4/NCV_1_ref_g4RNAFOLD_allG4COSMICfasta_var_features.txt")
alt_RNAFOLD_ncv<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/RNAfold/NCV_1_alt_g4/NCV_1_alt_g4RNAFOLD_allG4COSMICfasta_var_features.txt")

RNAFOLD_NCV<-left_join(ref_RNAFOLD_ncv[,c("MFE","MFEadj","EFE","CDE","MEAFE","MEA","ED","MFEadj.GC","MFEadj.BP","MEAFEadj","EDadj","CFE","EFEadj","name")],alt_RNAFOLD_ncv[,c("MFE","MFEadj","EFE","CDE","MEAFE","MEA","ED","MFEadj.GC","MFEadj.BP","MEAFEadj","EDadj","CFE","EFEadj","name")],by="name")
RNAFOLD_NCV<-RNAFOLD_NCV%>% mutate(file="NCV")



NCV_total<-left_join(NCV,NCVseq,by=c("name","G4_ref","G4_alt","G4_id"))
rm(NCV)
NCV1<-NCV_total
NCV1<-NCV1%>% filter(V12==1 | V12==-1 | V16==1 | V16== -1)
NCV1<-setDT(NCV1)
NCV1[, G4_refrc := ifelse((V12 == -1 | V16 == -1) , sapply( G4_ref, function(x) as.character(
reverseComplement(DNAString(x)))), sapply(G4_ref,function(x) as.character(x)))]


NCV1[, G4_altrc:= ifelse(V12 == -1 | V16 == -1, sapply( G4_alt, function(x) as.character(
reverseComplement(DNAString(x)))), sapply(G4_alt,function(x) as.character(x) ))]

NCV1$nameimp<-NCV1$name
NCV1$name<-gsub("_G4_id","",NCV1$name)
NCV1<-separate(data = NCV1, col = name, into = c("left", "right"), sep = "_")

NCV1$start<-as.numeric(NCV1$right)-30
NCV1$end<-as.numeric(NCV1$right)+30
#NCV1$name<-paste0(NCV1$left,":",NCV1$start,"-",NCV1$end)
NCV1$name<-paste0(NCV1$left,":",NCV1$start,"-",NCV1$end,"||","loc:",NCV1$loc_start,"||",NCV1$REF,"&",NCV1$ALT)
NCV2<-NCV1

writeFastaref(NCV2,"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/NCV/NCV_1_ref_g4.fasta")
writeFastaalt(NCV2,"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/NCV/NCV_1_alt_g4.fasta")


save(coding2,NCV2,file="/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/full_codingNCV_G4_duplicatesremoved.Rdata")
#load("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/full_codingNCV_G4_duplicatesremoved.Rdata")




coding2<- unique(coding2, by = c('name','V9','V12','V13','V16'))
NCV2<- unique(NCV2, by = c('name','V9','V12','V13','V16'))


value<-function(ref_RNAFOLD,alt_RNAFOLD,g4hunter2){
ref_RNAFOLD<-setDT(ref_RNAFOLD)
alt_RNAFOLD<-setDT(alt_RNAFOLD)
ref_RNAFOLD<- unique(ref_RNAFOLD, by = c('name'))

alt_RNAFOLD<- unique(alt_RNAFOLD, by = c('name'))
RNAFOLD<-merge(ref_RNAFOLD[,c("ID","MFE","MFEadj","EFE","CDE","MEAFE","MEA","ED","MFEadj.GC","MFEadj.BP","MEAFEadj","EDadj","CFE","EFEadj","name")],alt_RNAFOLD[,c("MFE","MFEadj","EFE","CDE","MEAFE","MEA","ED","MFEadj.GC","MFEadj.BP","MEAFEadj","EDadj","CFE","EFEadj","name")],by="name")
head(RNAFOLD)

#
whichdiff<-function(df,col){
df_orig<-df
colref<-paste0(col,".x")
colalt<-paste0(col,".y")
colnamedf<<-c(colref,colalt)
print(colnamedf)
df<-as.data.frame(df[,colnamedf])
df[,"which"]<-ifelse(df[,colref] > df[,colalt],"alt_sml","ref_sml")
df[,"which"]<-ifelse(df[,colref] ==df[,colalt],"no change",df[,"which"])
print(table(df$which))
return(cbind(df_orig,df[!names(df) %in% colnamedf]))
}

df1<-whichdiff(as.data.frame(data_dfannot1),"MFE")
levels(data_dfannot1$which)<-c("Further stabilized","no change","Destabilized")
df1$which<-as.factor(df1$which)
data_dfannot1$whichMFE<-(df1$which)

df1<-whichdiff(as.data.frame(data_dfannot1),"ED")
levels(df1$which)<-c("less diversity","no change","more diversity")
df1$which<-as.factor(df1$which)
data_dfannot1$whichED<-(df1$which)


#whichdiff(as.data.frame(RNAFOLD),"MFEadj")
#whichdiff(as.data.frame(RNAFOLD),"EFE")
#whichdiff(as.data.frame(RNAFOLD),"CDE")
#whichdiff(as.data.frame(RNAFOLD),"MEAFE")
#whichdiff(as.data.frame(RNAFOLD),"MEA")
#ED is positive so result will be opposite
#whichdiff(as.data.frame(RNAFOLD),"ED")
#whichdiff(as.data.frame(RNAFOLD),"CFE")
RNAFOLD<-df1
RNAFOLD$which<-as.factor(RNAFOLD$which)

levels(data_dfannot1$whichMFE)<-c("Further stabilized","no change","Destabilized")
#for ED
#levels(data_dfannot1$whichED)<-c("less diversity","no change","more diversity")


RNAFOLD<-separate(
  RNAFOLD,
  col=name,
  into=c("region","loc","SNP"),
  sep = "\\|\\|",
  remove = FALSE)
head(RNAFOLD)

RNAFOLD$SNP<-gsub("&","->",RNAFOLD$SNP)

merged<-merge(g4hunter2,RNAFOLD,by="name")
table(paste0(merged$V12,"::",merged$V16),merged$which)




merged$SNP <- as.character(merged$SNP)

test<- merged%>% 
       mutate(SNP_strand = case_when((V10 == 1 | V16 == 1) ~ merged$SNP,
                                         (V10 == -1 | V16 == -1) ~ chartr('ATGC', 'TACG',merged$SNP) ,
                                         (V10 == 1 & V16 == -1) ~ merged$SNP,
                                         (V10 == -1 & V16 == 1) ~ chartr('ATGC', 'TACG',merged$SNP))) 



test<-test%>% 
       mutate(G_quadstrand = case_when((V10 == 1 | V16 == 1) ~ "+",
                                         (V10 == -1 | V16 == -1) ~ "-" ,
                                         (V10 == 1 & V16 == -1) ~ "+",
                                         (V10 == -1 & V16 == 1) ~ "-")) 

return(test)
}


ref_RNAFOLD_coding<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/RNAfold/coding_1_ref_g4/coding_1_ref_g4RNAFOLD_allG4COSMICfasta_var_features.txt")
alt_RNAFOLD_coding<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/RNAfold/coding_1_alt_g4/coding_1_alt_g4RNAFOLD_allG4COSMICfasta_var_features.txt")
ref_RNAFOLD_NCV<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/RNAfold/NCV_1_ref_g4/NCV_1_ref_g4RNAFOLD_allG4COSMICfasta_var_features.txt")
alt_RNAFOLD_NCV<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/RNAfold/NCV_1_alt_g4/NCV_1_alt_g4RNAFOLD_allG4COSMICfasta_var_features.txt")



test_coding<-value(ref_RNAFOLD_coding,alt_RNAFOLD_coding,coding2)
test_NCV<-value(ref_RNAFOLD_NCV,alt_RNAFOLD_NCV,NCV2)

data<-list(CODING=test_coding,NCV=test_NCV)


#-------------------------------------
 NCV_total<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/CosmicNCV.tsv")

library(stringr)

NCV<-NCV_total %>%
  # Use hg38 coordinates and only consider frameshift (Indel) or substitutions
  filter(
    GRCh == 38 ) %>%
  # Extract chromosome coordinates from genome position' string
  mutate(
    chr   = str_c('chr', str_extract(`genome position`, '.*(?=:)')),
    start = str_extract(`genome position`, '(?<=:).*(?=-)') %>% as.integer(),
    end   = str_extract(`genome position`, '(?<=-).*') %>% as.integer()
  )

# Annotate with Gene information from UCSC using genome coordinates
NCV_annotated <-
 NCV  %>%
  select(Sampleid= `Sample name`,
	ref= WT_SEQ,
	alt= MUT_SEQ,
    primary_site      = `Primary site`,
    primary_histology = `Primary histology`,
    #zygosity          = `Mutation zygosity`,
    #somatic           = `Mutation somatic status`,
    mutation_id       = `GENOMIC_MUTATION_ID`,
   # mutation_desc     = `Mutation Description`,
    #mutation_cds      = `Mutation CDS`,
    #mutation_aa       = `Mutation AA`,
    fathmm_score= FATHMM_MKL_NON_CODING_SCORE,
    chr,
    #strand            = `Mutation strand`,
    start,
    end, 
  Mutation.genome.position= `genome position`) %>%
 distinct %>%
  mutate(
    cancer_type = case_when(
      str_detect(.$primary_histology, 'melanoma')                      ~ 'Malignant melanoma',
      .$primary_site == 'large_intestine'                              ~ 'Colorectal',
      .$primary_site == 'endometrium'                                  ~ 'Endometrial',
      .$primary_site == 'lung'                                         ~ 'Lung',
      .$primary_site == 'liver'                                        ~ 'Liver',
      .$primary_site == 'skin'                                         ~ 'Non-melanoma skin',
      .$primary_site == 'breast'                                       ~ 'Breast',
      .$primary_site == 'stomach'                                      ~ 'Stomach',
      .$primary_site %in% c('upper_aerodigestive_tract', 'oesophagus') ~ 'Upper aerodigestive',
      .$primary_site == 'haematopoietic_and_lymphoid_tissue'           ~ 'Blood',
      .$primary_site == 'prostate'                                     ~ 'Prostate',
      .$primary_site == 'pancreas'                                     ~ 'Pancreatic',
      .$primary_site == 'urinary_tract'                                ~ 'Bladder',
      .$primary_site == 'kidney'                                       ~ 'Kidney',
      .$primary_histology == 'glioma'                                  ~ 'Glioma',
      .$primary_site == 'ovary'                                        ~ 'Ovarian',
      .$primary_site == 'cervix'                                       ~ 'Cervical',
      .$primary_site == 'thyroid'                                      ~ 'Thyroid',
      .$primary_site == 'bone'                                         ~ 'Bone',
      .$primary_site == 'thymus'                                       ~ 'Thymus',
      .$primary_site == 'testis'                                       ~ 'Testicular',
      .$primary_site == 'small_intestine'                              ~ 'small Intestine',
      .$primary_site == 'central_nervous_system'                       ~ 'CNS',
      .$primary_site == 'biliary_tract'                       	     ~ 'Bile duct',
      TRUE ~ 'Other'))


#COSMIC_annotate
#here we make NCV_annotated and coding(COSMIC_annotated) same columns
NCV_annotated$mutation_cds<-"NH"
NCV_annotated$mutation_aa<-"NH"
NCV_annotated$HGVSG<-"NH"

NCV_annotated$strand<-"NH"
NCV_annotated$mutation_desc<-"NH"

#---------
#analyse for specific cancer_type

all<-rbind(COSMIC_annotated,NCV_annotated)
all<-all%>% filter(ref!="" | alt!="")

# load("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/RESULTS_final/COSMIC_geo_expr_g4.RData")
load("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/RESULTS_final/COSMIC_geo_expr_g4_annotated.RData")

G4_var_GSM<-G4_expr
colnames(G4_var_GSM)[1:3]<-c("chr","start","end")

G4_var_GSM<-G4_var_GSM %>% mutate(mut=gsub("^.*[||]", "", id))%>% separate(mut, into = c("ref","alt"))

sth<-left_join(G4_var_GSM,all,by=c("chr","start","end","ref","alt"))

temp<-sth%>% filter(mutation_desc!="Substitution - Missense")%>%filter(cancer_type=="Glioma")%>%
select(typemut,id,Sampleid,SNP_strand,cancer_type)%>%
 distinct()%>%group_by(cancer_type,SNP_strand)%>%
summarise(n=n())%>% as.data.frame()

temp<-sth%>%filter(cancer_type=="Glioma")%>%
select(typemut,id,Sampleid,SNP_strand,cancer_type,which1)%>% 
 distinct()%>%group_by(SNP_strand,id,which1)%>%
summarise(n=n(),totals=n_distinct(Sampleid),totalid=n_distinct(id))%>% arrange(desc(totals))%>%as.data.frame()

temp<- sth%>% filter(!is.na(Sampleid))%>%select(id,Sampleid,cancer_type)%>% distinct()%>%
select(cancer_type,Sampleid)%>% distinct()

glioma_TCGA_sample_id<-sth %>%filter(cancer_type=="Glioma")%>%
select(typemut,id,Sampleid,SNP_strand,cancer_type)%>%
  filter(stringr::str_detect(Sampleid, 'TCGA') )%>% group_by(Sampleid)%>% summarise(n=n())%>% arrange(desc(n))%>% as.data.frame()%>% head()


#-------------------------------------------

NCV_total<-separate(NCV_total,col="genome position",sep=":",into=c("chr","pos"),remove=FALSE)


NCV_total$chr<-paste0("chr",NCV_total$chr)

NCV_total$position<-str_split(NCV_total$pos, "\\-", simplify=T)[,1]


ncvbed<-NCV_total%>% dplyr::select(chr,position,position,`Primary site`,`Sample name`,`Primary histology`)%>%
summarise(chr=chr,start=position,end=position,name=paste0(`Sample name`,"~",`Primary site`,"~",`Primary histology`),score=0,strand="+")

coding_total$Mutation<-str_right(coding_total$HGVSG,3)

coding_total<-separate(coding_total,col="Mutation",sep=">",into=c("REF","ALT"),remove=TRUE)

codingbed<-coding_total%>% dplyr::select(chr,position,position,`Primary site`,`Sample name`,`Primary histology`,`Mutation strand`)%>%
summarise(chr=chr,start=position,end=position,name=paste0(`Sample name`,"~",`Primary site`,"~",`Primary histology`),score=0,strand=`Mutation strand`)


temp<-rbind(codingbed,ncvbed)
temp<-temp%>% filter(!(is.na(start) | is.na(end)))%>%summarise(chr=chr,start=(as.numeric(start)-1),end=end,name=name,score="0",strand=strand)

fwrite(temp,file="/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/TOTAL_Mutation.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

NCV<-NCV_total[,c("Sample name","Primary site","Primary histology",
"GENOMIC_MUTATION_ID","LEGACY_MUTATION_ID","genome position",
"WT_SEQ","MUT_SEQ","HGVSG")]
NCV<-separate(NCV,col="genome_position",sep=":",into=c("chr","pos"),remove=FALSE)
colnames(NCV)<-c("Sample_ID","Primary_site","Primary_histology",
"GENOMIC_MUTATION_ID","LEGACY_MUTATION_ID","genome_position",
"WT_SEQ","MUT_SEQ","HGVSG")

NCV<-NCV%>%mutate(site_histology=paste0(Primary_site,"~",Primary_histology))

NCV$chr<-paste0("chr",NCV$chr)

NCV$position<-str_split(NCV$pos, "\\-", simplify=T)[,1]

NCV$Primary_site<-NULL
NCV$Primary_histology<-NULL

NCV_merged<-merge(data[["NCV"]],NCV,by.x=c("REF","ALT","right"),by.y=c("WT_SEQ","MUT_SEQ","position"))
#------------------------------------------




 coding_total<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/MutantExport.tsv")


coding<-coding_total %>%
  # Use hg38 coordinates and only consider frameshift (Indel) or substitutions
  filter(
    GRCh == 38 ) %>%
  # Extract chromosome coordinates from 'Mutation genome position' string
  mutate(
    chr   = str_c('chr', str_extract(`Mutation genome position`, '.*(?=:)')),
    start = str_extract(`Mutation genome position`, '(?<=:).*(?=-)') %>% as.integer(),
    end   = str_extract(`Mutation genome position`, '(?<=-).*') %>% as.integer()
  )


    # Remove duplicate positions
	coding$`Gene name` <- gsub("_ENST\\d+", "", coding$`Gene name`)
	nrow(coding)
	coding<-as.data.frame(coding)
   	coding<- coding[!duplicated(coding[c("Gene name",  "mutation genome position")]), ]
	nrow(coding)
str_right <- function(string, n) {
  substr(string, nchar(string) - (n - 1), nchar(string))
}


# Annotate with Gene information from UCSC using genome coordinates
COSMIC_annotated <-
 coding  %>%
  select(Sampleid= `Sample name`,
    primary_site      = `Primary site`,
    primary_histology = `Primary histology`,
    #zygosity          = `Mutation zygosity`,
    #somatic           = `Mutation somatic status`,
    mutation_id       = `MUTATION_ID`,
    mutation_desc     = `Mutation Description`,
    mutation_cds      = `Mutation CDS`,
    mutation_aa       = `Mutation AA`,
	fathmm_score=`FATHMM score`,
    chr,
    strand            = `Mutation strand`,
    start,
    end, 
    HGVSG,
  Mutation.genome.position= `Mutation genome position`) %>%
 distinct %>%
  mutate(Mutation=str_right(HGVSG,3),
    cancer_type = case_when(
      str_detect(.$primary_histology, 'melanoma')                      ~ 'Malignant melanoma',
      .$primary_site == 'large_intestine'                              ~ 'Colorectal',
      .$primary_site == 'endometrium'                                  ~ 'Endometrial',
      .$primary_site == 'lung'                                         ~ 'Lung',
      .$primary_site == 'liver'                                        ~ 'Liver',
      .$primary_site == 'skin'                                         ~ 'Non-melanoma skin',
      .$primary_site == 'breast'                                       ~ 'Breast',
      .$primary_site == 'stomach'                                      ~ 'Stomach',
      .$primary_site %in% c('upper_aerodigestive_tract', 'oesophagus') ~ 'Upper aerodigestive',
      .$primary_site == 'haematopoietic_and_lymphoid_tissue'           ~ 'Blood',
      .$primary_site == 'prostate'                                     ~ 'Prostate',
      .$primary_site == 'pancreas'                                     ~ 'Pancreatic',
      .$primary_site == 'urinary_tract'                                ~ 'Bladder',
      .$primary_site == 'kidney'                                       ~ 'Kidney',
      .$primary_histology == 'glioma'                                  ~ 'Glioma',
      .$primary_site == 'ovary'                                        ~ 'Ovarian',
      .$primary_site == 'cervix'                                       ~ 'Cervical',
      .$primary_site == 'thyroid'                                      ~ 'Thyroid',
      .$primary_site == 'bone'                                         ~ 'Bone',
      .$primary_site == 'thymus'                                       ~ 'Thymus',
      .$primary_site == 'testis'                                       ~ 'Testicular',
      .$primary_site == 'small_intestine'                              ~ 'small Intestine',
      .$primary_site == 'central_nervous_system'                       ~ 'CNS',
      .$primary_site == 'biliary_tract'                       	     ~ 'Bile duct',
      TRUE ~ 'Other'))

#coding$Mutation<-str_right(coding$HGVSG,3)

COSMIC_annotated$check<-substring(COSMIC_annotated$Mutation, 2, 2)==">"
nrow(COSMIC_annotated)
COSMIC_annotated<-COSMIC_annotated %>% filter(check==TRUE)
nrow(COSMIC_annotated)
COSMIC_annotated <-separate(COSMIC_annotated,col="Mutation",sep=">",into=c("ref","alt"),remove=TRUE)
COSMIC_annotated$check<-NULL

#----------------------------------------------------


coding_total<-separate(coding_total,col="Mutation genome position",sep=":",into=c("chr","pos"),remove=FALSE)


coding_total$chr<-paste0("chr",coding_total$chr)

coding_total$position<-str_split(coding_total$pos, "\\-", simplify=T)[,1]




coding<-coding_total[,c("Gene name","Sample name","HGNC ID","Primary site","Primary histology",
"GENOMIC_MUTATION_ID","LEGACY_MUTATION_ID","Mutation genome position","Mutation CDS", "Mutation AA",
"Mutation Description",
"Mutation strand","FATHMM prediction","FATHMM score","Tumour origin",
"HGVSG")]

colnames(coding)<-c("gene","Sample_ID","HGNC_id","Primary_site","Primary_histology",
"GENOMIC_MUTATION_ID","LEGACY_MUTATION_ID","mutation_genome_position",
"Mutation_CDS","Mutation_AA",
"Mutation_Description",
"Mutation.strand","FATHMM_prediction","FATHMM_score","Tumour_origin","HGVSG")

    names(coding) <- tolower(gsub("\\.", "_", names(coding)))

    # Simplify COSMIC sample/site names
    coding$sample_id <- toupper(gsub("[-. ]", "", coding$sample_id))
    coding$primary_site <- toupper(gsub("[-. ]", "", coding$primary_site))
    coding<- coding[coding$mutation_genome_position != "", ]


 # Keep only SNVs
toMatch<-c("Substitution","Unknown")

library(data.table)
coding[coding$mutation_description %like% (paste(toMatch,collapse="|")),]

    coding1<- coding[grep(pattern='Substitution | Unknown', coding$mutation_description), ]
matches <- unique (grep(paste(toMatch,collapse="|"), 
                        coding$mutation_description, value=FALSE))
	nrow(coding)
    # Separate chromosomes and positions
    positions <- data.frame(do.call(rbind, strsplit(as.vector(
        coding$mutation_genome_position), split = ":|-")))
    names(positions) <- c("chr", "start", "end")
    coding<- cbind(coding, positions)
    coding$start <- as.numeric(as.character(coding$start))
    coding$end <- as.numeric(as.character(coding$end))
    coding$mutation_genome_position <- NULL

coding<-separate(coding,col="genome_position",sep=":",into=c("chr","pos"),remove=FALSE)

coding<-coding%>%mutate(site_histology=paste0(Primary_site,"_",Primary_histology))



coding$chr<-paste0("chr",coding$chr)

coding$position<-str_split(coding$pos, "\\-", simplify=T)[,1]
rm(coding_total)


coding$site_histology<-gsub("_","~",coding$site_histology)

coding$col<-paste0(coding$Sample_ID,"~",coding$site_histology)





coding$Primary_site<-NULL
coding$Primary_histology<-NULL
coding$pos<-NULL
coding$Mutation<-str_right(coding$HGVSG,3)
coding<-separate(coding,col="Mutation",sep=">",into=c("REF","ALT"),remove=TRUE)

coding_merged<-merge(data[["CODING"]],coding,
by.x=c("REF","ALT","right"),by.y=c("REF","ALT","position"))



#starting point for any analysis
#--------------------------------------
save(coding_merged,NCV_merged,file="/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/CODING_NCV_merged.RData")



map(data,.x%>% 
#load("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/CODING_NCV_merged.RData")

data<-list(coding=coding_merged,NCV=NCV_merged)


data[[1]]<-data[[1]]%>% select(-c(gene))%>% distinct()  
#data[[2]]<-data[[2]]%>% select(-c(gene))%>% distinct()  

data[[2]]$pos<-NULL
#c(name,G4_id,G4_refrc,G4_altrc,V9,V12,V13,V16,loc_start,region,loc_start,SNP_strand,
data<-map(data,~.x%>% select(-c(REF,ALT,right,G4_id,G4_ref,G4_alt,V10,V11,V14,V15,seq_g4_full,nameimp,ID,loc,SNP,
chr,Sample_ID,site_histology,genome_position,GENOMIC_MUTATION_ID,LEGACY_MUTATION_ID,HGVSG))%>% distinct())

##here we have to split some annotations, like sample_id, and site histology. they arent relevant for this because we are just accounting
#for one SNP each position based on the alt allele. if multiple substitutions in different, we keep
for sample level study we come back to above data.



#data<-map(data,~setkey(.x, left, start, end))
data_df<-rbindlist(data,idcol="source")
data_df<-data_df%>% arrange(source)%>%group_by(left,loc_start,SNP_strand)%>%  
mutate(source=paste0(unique(source),collapse="|"))%>%distinct()%>%data.table()

#data_df%>% summarise(left,loc_start,end=loc_start,name,SNP_strand,G_quadstrand)%>% fwrite(.,file="cosmic_G4_variants.bed",sep="\t",col.names=FALSE)

data_df%>% group_by(left,loc_start,SNP_strand)%>%  count()%>% arrange(desc(n))


#-------------------------------------------------
#check if mutation lies inside G4 or flanking
temp_coding_r<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/coding/coding_1_ref_g4_quadparser.tsv",sep="\t")
temp_coding_a<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/coding/coding_1_alt_g4_quadparser.tsv",sep="\t")

temp_NCV_r<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/NCV/NCV_1_ref_g4_quadparser.tsv",sep="\t")
temp_NCV_a<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/NCV/NCV_1_alt_g4_quadparser.tsv",sep="\t")

temp_var_r<-rbind(temp_coding_r,temp_NCV_r)
temp_var_a<-rbind(temp_coding_a,temp_NCV_a)
colnames(temp_var_r)[c(2,3,4,5)]<-c("g4_start","g4_stop","g4seq","g4_type")
colnames(temp_var_a)[c(2,3,4,5)]<-c("g4_start","g4_stop","g4seq","g4_type")
temp_var_a%>% group_by(V1)%>% mutate(n=n())%>% arrange(desc(n))%>% head()
temp_var_a<-temp_var_a[!duplicated(temp_var_a), ]
temp_var_r<-temp_var_r[!duplicated(temp_var_r),]
rm(temp_coding_r,temp_coding_a,temp_NCV_a,temp_NCV_r)

data1<-left_join(data_df,temp_var_r[,c("V1","g4_start","g4_stop","g4seq","g4_type")],by=c("name"="V1"))
data1na<-data1%>% filter(is.na(g4_start))%>% select(-c(g4_start,g4_stop,g4seq,g4_type))
data1<-data1%>% filter(!is.na(g4_start))

data1na<-left_join(data1na,temp_var_a[,c("V1","g4_start","g4_stop","g4seq","g4_type")],by=c("name"="V1"))

data_total<-rbind(data1,data1na)
rm(data1  ,data1na)

data_total_1<-data_total%>% as.data.frame()%>%group_by(name,G4_refrc,G4_altrc,V12,V16,start,G_quadstrand,SNP_strand,g4_type)%>% 
mutate(numberofg4=n(),which1=paste0(unique(which),collapse="|"),g4seq1=paste0(unique(g4seq),collapse="|"),g4_start=min(as.numeric(g4_start)),g4_stop=max(as.numeric(g4_stop)))%>%
arrange(desc(numberofg4)) %>%
 mutate(g4seq=g4seq1,g4seq1=NULL,which=which1,which=NULL)%>%distinct()%>% as.data.frame()


data_total_2<-data_total_1%>% 
       mutate(affectg4 = case_when(((g4_start < 30) & (g4_stop<30)) ~"mut_downstream",
                                        ((g4_start >30) & (g4_stop>31)) ~ "mut_upstream", 
					((g4_start <31) & (g4_stop>29)) ~ "affect_g4"))
						

data_total_1<-data_total_2%>% 
       mutate(rel_loc= case_when((affectg4 =="mut_downstream")~(30-g4_stop)/30,
                                        (affectg4 == "mut_upstream")~(30-g4_start)/30,
					(affectg4=="affect_g4")~ (30-g4_start)/(g4_stop-g4_start)))
save(data_total_1,file="/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/full_codingNCV_G4_duplicatesremoved.Rdata")
#load("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/full_codingNCV_G4_duplicatesremoved.Rdata")
data_total_2<-data_total_1%>% filter(affectg4=="affect_g4")%>% distinct()

#data_total_1%>% group_by(name,G4_refrc,G4_altrc,start,g4seq,rel_loc)%>% mutate(affectg41=paste0(unique(affectg4),collapse="|"),g4seq11=paste0(unique(g4seq),collapse="|"))%>% arrange(desc(numberofg4))%>% head()%>% as.data.frame() 
#data_total_1<-data_total_1[!duplicated(data_total_1$name),]

#data_total_1%>% group_by(name,G4_refrc,G4_altrc,start,end)%>% mutate(n=n())%>% arrange(desc(n))%>%head()

temp<-split(data_total_1,data_total_1$affectg4)

temp<-split(data_total_2,data_total_2$G_quadstrand)

map(temp,~hist(.x$rel_loc))
#data_total_1$rel_loc<-(data_total_1$loc_start-30+data_total_1$g4_start-data_total_1$loc_start)/nchar(data_total_1$g4seq)
#data_dfannot_1$rel_loc<-abs(data_dfannot_1$loc_start+30+data_dfannot_1$g4_start-data_dfannot_1$loc_start)/data_dfannot_1$g4_length
--------------------------------------------------------------
#ANALYSIScb

data_total_2$diff=data_total_2$MFE.y-data_total_2$MFE.x

data_total_2$gquad<-unlist(lapply(strsplit(as.character(data_total_2$g4_type), ":"), '[[', 1))

x<-data_total_2$diff[data_total_2$SNP_strand=="A->G"]
y<-data_total_2$diff[data_total_2$SNP_strand=="T->G"]




 dt<-data_total_2%>% filter(gquad ==  "4" | gquad== "5" |  gquad== "6")

 dt$gquad<-as.factor(dt$gquad)
 dt$which1<-as.factor(dt$which1)
dt$SNP_strand<-as.factor(dt$SNP_strand)



data_total_2%>% filter(affectg4=="affect_g4")%>%mutate(real_mut=gsub("^.*[||]", "", name))%>% select(real_mut,G_quadstrand)%>% group_by(real_mut,G_quadstrand)%>% 
summarise(n=n()) %>% as.data.frame()%>% cb(.,sep=",")

mat<-map(temp,~xtabs(~paste0(SNP_strand,"::")+paste0(V12,":",V16),.x))
mat<-map(temp,~xtabs(~left+paste0(SNP_strand,"::"),.x))
temp<-map(temp,~setDT(.x))
count_temp<-map(temp,~.x[, .(G4hunter_ref = mean(V9), G4hunter_alt = mean(V13),count=.N,sd=sd(V9),sd2=sd(V13)), by = .(V12,V16)])

map(temp,~.x[ , .N, by = .(V8, SNP_strand) ]%>% dcast(SNP_strand~V8,value.var="N"))#%>% cb(.)

library(gmodels)

map(data,~CrossTable(.x$SNP_strand,paste0(.x$V12,":",.x$V16), expected = TRUE))

tempmat<-map(temp,~xtabs(~paste0(SNP,"::")+paste0(V12,":",V16),.x))
map(data,~xtabs(~V12,V16,.x))

c<-map(mat,~chisq.test(.x[,c(1,7)]))
#also we get one g4 location per mutation
data_total_2%>% mutate(chr=gsub(":.*$", "", name),start=start+g4_start,end=start+g4_start+g4_stop,score=".")%>%
select(chr,start,end,name,score,G_quadstrand)%>%
fwrite(.,"G4_aroundSNVcosmic_all.bed",sep="\t",col.names=FALSE,quote=FALSE)


#-------------------------------------

cat G4_aroundSNVcosmic_all.bed | bedtools slop -b 2000 -i - -g hg38.sizes |  sort-bed - |
bedtools merge -s -i - -c 4,6 -o collapse,distinct -delim "+" | 
bedops --chop 100 - > G4_aroundSNVcosmic_flank_window100.bed



sort -k1,1n -k2,2n G4_aroundSNVcosmic_all.bed  > G4_aroundSNVcosmic_all.sorted.bed

#bedtools merge -s  -i G4_aroundSNVcosmic_all.sorted.bed -c 4,6 -o collapse,distinct -delim "+"  > G4_aroundSNVcosmic_all.sortedmerged.bed

bedtools makewindows -g ./hg38.sizes -w 1 | bedtools map -a - -b /bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/G4_aroundSNVcosmic_all.sortedmerged.bed -c 4 -o count > G4_aroundSNVcosmic_all.sortedmerged.counts_SNV.bed
 
#--------------------------------------------------------------

#SNV affecting the protein coding genes
temp_df<-data_dfannot%>% filter(V19=="protein-coding")
 temp_df[, sum := .N, by = V8][, prop := .N, by = list(SNP_strand, which,V8) ][, prop := prop/sum][, sum := NULL]%>% dcast(SNP_strand~which,value.var="N")
round(prop.table(table( temp_df$V8, temp_df$SNP_strand),1),2)
temp_df[ , .N, by = .(V8, SNP_strand) ]%>% dcast(SNP_strand~V8,value.var="N")#%>% cb(.)

df1<-data_dfannot

mat<-xtabs(~SNP+which,coding_merged)
mat<-xtabs(~paste0(V12,":",V16,":",which)+paste0(SNP,"::"),coding_merged)
tab2<-chisq.test(t(mat[c(1,2,3,18,19,20),]), simulate.p.value = T)

chisq.posthoc.test(t(mat[c(1,2,3,18,19,20),]),method = "bonferroni")
#%>% write.table(., "clipboard-16384", sep=",", row.names=FALSE,quote=FALSE)



#---


maf_coding<-coding_merged%>% 
dplyr::select(chr, right,right,REF,ALT,Sample_ID)%>% 
summarise(chr=chr,start=right,end=right,REF,ALT,Sample_ID)
maf_NCV<-NCV_merged%>% 
dplyr::select(chr, right,right,REF,ALT,Sample_ID)%>%
 summarise(chr=chr,start=right,end=right,REF,ALT,Sample_ID)

maf_merged<-list(coding=maf_coding,NCV=maf_NCV)
save(maf_merged,file="/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/maf_merged.RData")

#-------

readMAF <- function(MAF,genome,...) {
    maf = read.delim(MAF,comment.char='#',...)
    
    rowRanges = VRanges(seqnames=Rle(maf$Chromosome),
        ranges=IRanges(start=maf$Start_position,end=maf$End_position),
        ref=as.character(maf$Reference_Allele),
        alt=as.character(maf$Tumor_Seq_Allele1),
        sampleNames=maf$Tumor_Sample_Barcode)
    return(rowRanges)
}



mat<-map(temp,~xtabs(~paste0(SNP_strand,"::")+paste0(V12,":",V16),.x))
library(gmodels)

map(data,~CrossTable(.x$SNP_strand,paste0(.x$V12,":",.x$V16), expected = TRUE))

tempmat<-map(temp,~xtabs(~paste0(SNP,"::")+paste0(V12,":",V16),.x))
map(data,~xtabs(~V12,V16,.x))

c<-map(mat,~chisq.test(.x[,c(1,7)]))


mat<-xtabs(~SNP+which,coding_merged)
mat<-xtabs(~paste0(V12,":",V16,":",which)+paste0(SNP,"::"),coding_merged)
tab2<-chisq.test(t(mat[c(1,2,3,18,19,20),]), simulate.p.value = T)

chisq.posthoc.test(t(mat[c(1,2,3,18,19,20),]),method = "bonferroni")
#%>% write.table(., "clipboard-16384", sep=",", row.names=FALSE,quote=FALSE)


select(MFE.x,MFE.y,SNP,loc,name,V12,V16,which)%>%

a<-coding_merged%>% dplyr::select(SNP,loc,region,V12,V16,REF,ALT)%>% distinct()%>%as.data.frame()%>%
  summarise(chrom="chromosome",id=gsub("chr(.+):.*", "\\1", region),location=gsub("loc:","",loc),ref=REF,alt=ALT,strand="1")%>%arrange(id,location)%>% 
  write.table(., "clipboard-16384", sep="\t", row.names=FALSE,quote=FALSE,col.names=FALSE)



library(data.table) # v1.9.7 (devel version)

setDT(coding_merged) # convert your dataframe into a data.table

# save files
  a[, fwrite(.SD, paste0("output", var1,".txt")), by = id]

a[, fwrite(copy(.SD)[, id := id] paste0("output", id,".txt"),sep="\t",col.names=FALSE), by = id]

a1<-split(a,a$id)
lapply(names(a1), function(x){fwrite(a1[[x]], file = paste("output", x, sep = ""), sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)})





#----------------------------------------

data[[1]]<-data[[1]]%>% select(-c(gene))%>% distinct()  
#data[[2]]<-data[[2]]%>% select(-c(gene))%>% distinct()  

data[[2]]$pos<-NULL
#c(name,G4_id,G4_refrc,G4_altrc,V9,V12,V13,V16,loc_start,region,loc_start,SNP_strand,
data<-map(data,~.x%>% select(-c(REF,ALT,right,G4_id,G4_ref,G4_alt,V10,V11,V14,V15,seq_g4_full,nameimp,ID,loc,SNP,
chr,Sample_ID,site_histology,genome_position,GENOMIC_MUTATION_ID,LEGACY_MUTATION_ID,HGVSG))%>% distinct())

##here we have to split some annotations, like sample_id, and site histology. they arent relevant for this because we are just accounting
#for one SNP each position based on the alt allele. if multiple substitutions in different, we keep
for sample level study we come back to above data.



data<-map(data,~setkey(.x, left, start, end))

annotation<-setkey(annotation, V1, V2, V3)
#annot_list<-map(annot,annot$ANNOTATION)
data_df<-rbindlist(data,idcol="source")
data_df<-data_df%>% arrange(source)%>%group_by(left,loc_start,SNP_strand)%>%  
mutate(source=paste0(unique(source),collapse="|"))%>%distinct()%>%data.table()

#data_total_2%>% summarise(left,loc_start,end=loc_start,name,SNP_strand,G_quadstrand)%>% fwrite(.,file="cosmic_G4_variants.bed",sep="\t",col.names=FALSE)

data_df%>% group_by(left,loc_start,SNP_strand)%>%  count()%>% arrange(desc(n))
#-----------------------------------------
annotatePeaks.pl ./cosmic_G4_variants.bed hg38 -annStats cosmic_variants_stats.tsv > cosmic_G4_variants_annotated.bed
 ./annotatePeaks.pl /bio/home/goonerrn/g_quadruplex/annotation_database/SNP/ClinVar/clinvar_20200203.bed hg38 -annStats /bio/home/goonerrn/g_quadruplex/annotation_database/SNP/ClinVar/Clinvar_stats.tsv > /bio/home/goonerrn/g_quadruplex/annotation_database/SNP/ClinVar/clinvar_variants_annotated.bed 
./annotatePeaks.pl /bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/G4_aroundSNVCLINVAR_all.bed hg38 -annStats /bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/Clinvar_stats.tsv > /bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/G4_aroundSNVCLINVAR_all_annotated.bed 

data_dfannot<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/cosmic_G4_variants_annotated.bed",skip=1)

##annotation by annotatr 
#-------------------------------------------
f<-load("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/COSMIC_mutations_G4_aroundannotated.RData")
cosmic_annotated_SNV[,c("strand","seqnames","start","end","width","change","annot.start","annot.id","annot.tx_id", "annot.end","annot.width")]<-list(NULL)
head(cosmic_annotated_SNV)
df<-cosmic_annotated_SNV%>%mutate(annot.p=case_when((annot.type=="hg38_lncrna_gencode")~"lncRNA GENCODE",
(annot.type=="hg38_enhancers_fantom")~"Enhancers",
(annot.type=="hg38_genes_1to5kb")~"Nearby1to5kb",
                                  (annot.type=="hg38_cpg_islands")~"CpG Islands",
               (annot.type=="hg38_genes_promoters")~"PROMOTERS",
(annot.type=="hg38_genes_introns")~"INTRON",
(annot.type=="hg38_genes_cds")~"CDS",
(annot.type=="hg38_genes_exons")~"EXON",
(annot.type=="hg38_genes_5UTRs")~"5' UTR",
(annot.type=="hg38_genes_3UTRs")~"3' UTR",
(annot.type=="hg38_genes_exonintronboundaries")~"EXON/INTRON Boundaries",
(annot.type=="hg38_genes_intronexonboundaries")~"INTRON/EXON Boundaries",
(annot.type=="hg38_genes_intergenic")~"Intergenic"))


df$annot.p<- factor(df$annot.p,
                    levels = c("CDS", "5' UTR", "3' UTR", "EXON", "INTRON", "EXON/INTRON Boundaries","INTRON/EXON Boundaries",
"PROMOTERS","Enhancers", "hg38_genes_firstexons", "CpG Islands","lncRNA GENCODE","Intergenic" ),ordered=TRUE)
head(df)

df1<-df%>% filter(annot.p %in% c("CDS", "5' UTR", "3' UTR", "EXON", "INTRON","PROMOTERS","Enhancers","CpG Islands","lncRNA GENCODE","Intergenic"))

df2<-df1 %>% group_by(id,annot.p) %>% fill(annot.gene_id,annot.symbol)%>% distinct()%>%as.data.frame()#%>%head()


#df<-separate(df,col="annot.id",into=c("type","geneid"),remove=FALSE,sep=":")


df3<-df2%>% group_by(id,annot.seqnames,annot.p,annot.strand)%>% 
mutate(symbol=paste(unique(ifelse(is.na(annot.symbol), "", annot.symbol)), collapse = "|"),
geneid=paste(unique(ifelse(is.na(annot.gene_id), "", annot.gene_id)), collapse = "|")) %>%
mutate(annot.symbol=NULL,annot.gene_id=NULL)%>% distinct()%>% as.data.frame()


temp<- df3%>%
   mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))%>% arrange(annot.p)

#temp<-df3 %>% arrange(annot.p)%>% as.data.frame()%>%group_by(id,geneid) %>%
#  dplyr::slice(1) %>%
#  ungroup()%>% as.data.frame()


a<-temp%>% arrange(symbol)%>%group_by(id)%>% filter(!(is.na(geneid) | annot.strand=="*"))%>% 
 summarise(n=n(),annot=paste(unique(annot.p),collapse="|"),gene.id=paste(unique(geneid),collapse="|"),genesymbol=paste(unique(symbol),collapse="|"))%>%
 arrange(desc(n))%>% as.data.frame()

all_genes_g4<-data_dfannot1%>% filter(!(annot.gene_id=="" | is.na(annot.gene_id))) %>% group_by(annot.gene_id,name)%>% summarise(n=n())%>% ungroup()%>% select(annot.gene_id)%>%  mutate(V2 = strsplit(as.character(annot.gene_id), "\\|")) %>% 
summarise(V3=strsplit(as.character(V2), " "))%>% select(V3)%>%
    unnest(V3) %>%  distinct()
head(temp)
(df3)%>%
filter(annot.p!="Intergenic" & annot.p!="lncRNA GENCODE" & !is.na(symbol))%>% head(20)


#-----------------------------------------------
 awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' /bio/data/gtf/Hs/TCGA/gencode.v22.annotation.gtf | gtf2bed  -  > GENCODE_annotation_hg38.bed  |
perl -lne 'print "@m" if @m=(/((?:gene_type|gene_name)\s+\S+)/g);' -  > GENCODE_annotation_hg38_parsetocbind.tsv


cut -f1-9 GENCODE_annotation_hg38.bed | paste --delimiters '\t'  -  GENCODE_annotation_hg38_parsetocbind.tsv > GENCODE_annotation_hg38f.bed

#-----------------------------------------------------------


data_dfannot1<-left_join(data_total_2,df2,by=c("name"="id"))

data_dfannot1%>% filter(diff > 10 | diff < -5)%>% filter(annot.p=="PROMOTER")%>% filter(which1=="Destabilized")%>% select(annot.gene_id) %>%head() 
data_dfannot1<-data_dfannot1%>% mutate(diff =(MFE.y-MFE.x))
data_dfannot1<-data_dfannot1%>% mutate(chrom=paste0("chr",gsub("chr(.+):.*", "\\1", region)))
data_dfannot1$mut<-sapply(strsplit(name, "||", fixed=TRUE), tail, 1)))


setkey(data_dfannot1,chrom,loc_start)


vcf$mut<-paste0(vcf$REF,"&",vcf$ALT)

temp<-merge(data_dfannot1, vcf, by.x=c('chrom', 'loc_start','mut'), by.y=c('CHROM', 'POS','mut'))
temp<-separate(temp,INFO,";",into=c("AC","AF"))

temp$AC<-(gsub("AC=","",temp$AC))
temp$AF<-(gsub("AF=","",temp$AF))

temp$AC<-sapply(strsplit(as.character(temp$AC),","), function(x) {
       x1 <- as.numeric(x)
       if(all(is.na(x1))) NA_real_ else max(x1, na.rm = TRUE)
  })


temp$AF<-sapply(strsplit(as.character(temp$AF),","), function(x) {
       x1 <- as.numeric(x)
       if(all(is.na(x1))) NA_real_ else max(x1, na.rm = TRUE)
  })


temp<-temp %>% 
  mutate(diff_per = ntile(diff, 100))

vcf$POS<-as.character(vcf$POS)
data_dfannot1$loc_start<-as.character(data_dfannot1$loc_start)


data_dfannot1[, as.list(summary(diff)), by = .(annot.p,which1)]
data_dfannot1%>% group_by(which1)%>%summarise(N.num=n(),
            mean.num = mean(diff),
            sd.num = sd(diff)) %>%
  mutate(se.num = sd.num / sqrt(N.num),
         lower.ci = 100*(mean.num - qt(1 - (0.05 / 2), N.num - 1) * se.num),
         upper.ci = 100*(mean.num + qt(1 - (0.05 / 2), N.num - 1) * se.num))


[, as.list(summary(diff)), by = .(annot.p,which1)]


table2<-data_dfannot1%>%filter(gquad %in% c(4,5)) %>% filter(!is.na(annot.p))%>%filter(which1!="no change")%>% select(annot.p,which1,gquad,name,diff)%>% distinct()%>%
  group_by(annot.p,which1,gquad) %>% 
  summarise(N.num=n(),
            mean.num = mean(diff),median=(median(diff)),q90 = quantile(diff, probs = 0.90), 
           std= sqrt(var(diff)),se.num = std/ sqrt(N.num),
            lower = mean(diff) - qnorm(.975)*std/sqrt(n()),
            upper = mean(diff) + qnorm(.975)*std/sqrt(n()))%>%as.data.frame()#%>% cb(.)

library(plyr) # function ddply()
mydf_m <- ddply(data_dfannot1, .(which1), transform, ecd = ecdf(abs(diff))(abs(diff)))

ggplot(mydf_m,aes(x = abs(diff), y = ecd)) + 
    geom_line(aes(group = which1, colour = which1))

data_dfannot1%>%  filter(diff > 10 | diff < -10)%>% filter(annot.p=="INTRON")%>%
 filter(which1=="Destabilized")%>% select(annot.gene_id) %>%head() 


data_dfannot1%>%  filter(which1=="Destabilized")%>%
filter(diff> quantile(diff,0.95))%>% head()


data_dfannot1$g4genestrand<-ifelse(data_dfannot1$annot.strand==data_dfannot1$G_quadstrand,"same","opposite")
data_dfannot1$g4genestrand<-ifelse(data_dfannot1$annot.strand=="*","notapplicable",data_dfannot1$g4genestrand)
transitions<-c("A->G","G->A","C->T","T->C")
transversions<-c("A->T","A->C","C->A","C->G","T->A","T->G","G->C","G->T")
data_dfannot1$typemut<-ifelse(data_dfannot1$SNP_strand %in% transitions,"Transition","Transversion")



#---------------------------------
#this might not work, check
# creating dataframe with mean diff per annotation
cty_mpg <- 
data_dfannot1  %>% filter(!is.na(annot.p))%>%filter(SNP_strand=="A->G" & which1!="no change")%>%
dplyr::group_by(.data = ., annot.p) %>%
  dplyr::summarise(.data = ., mean_diff = median(diff, na.rm = TRUE)) %>%
  dplyr::rename( make = annot.p) %>%
  dplyr::arrange( mean_diff) %>%
  dplyr::mutate( make = factor(x = make, levels = .$make)) %>%
  dplyr::mutate(
    .data = .,
    percent_rank = (trunc(rank(mean_diff)) / length(mean_diff)) * 100
  ) %>%distinct()%>%
  tibble::as_data_frame(x = .)


# plot
ggplot2::ggplot(data = cty_mpg, mapping = ggplot2::aes(x = make, y = mean_diff)) +
  ggplot2::geom_point(col = "tomato2", size = 3) + # Draw points
  ggplot2::geom_segment(
    mapping = ggplot2::aes(
      x = make,
      xend = make,
      y = min(mean_diff),
      yend = max(mean_diff)
    ),
    linetype = "dashed",
    size = 0.1
  ) + theme_bw()+# Draw dashed lines
  ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(trans = ~(trunc(rank(.)) / length(.)) * 100, name = "percentile")) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = expression("Effect of A->G SNV on"~Delta~MFE), 
       x = "Annotation",
	y = expression(Delta ~MFE),
    subtitle = "Dot plot (A->G SNV) ",
    caption = "source: A->G variant non 0"
  ) 
dev.off()

labs(title=expression("Effect of T->G SNV on"~Delta~MFE),
x = "Relative Location of SNV in a G4 region", 
       y = expression(Delta ~MFE)) + 
  theme_minimal() +

#------------------------------------------------------------------------

cty_mpg <- 
data_dfannot1%>% filter(!is.na(annot.p))%>% filter(SNP_strand=="T->G"  & which1!="no change")%>%
dplyr::group_by(.data = ., annot.p) %>%
  dplyr::summarise(.data = ., mean_diff = median(diff, na.rm = TRUE)) %>%
  dplyr::rename( make = annot.p) %>%
  dplyr::arrange( mean_diff) %>%
  dplyr::mutate( make = factor(x = make, levels = .$make)) %>%
  dplyr::mutate(
    .data = .,
    percent_rank = (trunc(rank(mean_diff)) / length(mean_diff)) * 100
  ) %>%
  tibble::as_data_frame(x = .)

# plot
ggplot2::ggplot(data = cty_mpg, mapping = ggplot2::aes(x = make, y = mean_diff,color=)) +
  ggplot2::geom_point(col = "tomato2", size = 3) + # Draw points
  ggplot2::geom_segment(
    mapping = ggplot2::aes(
      x = make,
      xend = make,
      y = min(mean_diff),
      yend = max(mean_diff)
    ),
    linetype = "dashed",
    size = 0.1
  )  + theme_bw()+## Draw dashed lines
  ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(trans = ~(trunc(rank(.)) / length(.)) * 100, name = "percentile")) +
  ggplot2::coord_flip() +
    ggplot2::labs(
    title = expression("Effect of T->G SNV on"~Delta~MFE), 
       x = "Annotation",
	y = expression(Delta ~MFE),
    subtitle = "Dot plot",
    caption = "source: T->G variant non 0 cosmic"
  )  
dev.off()

#-------------------------------------------------------------------------

#_------------------------------------------------------------------------

temp1<-data_dfannot1%>% filter(!is.na(annot.p) & which1!="no change")%>%
  filter(SNP_strand=="T->G")%>% select(diff,name,annot.p,which1)%>% distinct()
temp1<- aggregate(diff~annot.p, temp1, FUN = median)
#temp1$diff<-temp1$diff*-1
temp1<- data.frame(temp1[order(temp1$diff), ], rank=1:nrow(temp1))

g <- ggplot(temp1, aes(y = rank, x = diff))
g <- g + geom_point(size = 3,col = "tomato3")
g <- g + scale_y_continuous(name = "Annotation", labels = temp1$annot.p, breaks = temp1$rank,
                            sec.axis = dup_axis(name = element_blank(),
                                                breaks = seq(1, nrow(temp1), (nrow(temp1)-1)/4),
                                                labels = 25 * 0:4))

#g <- g + scale_x_continuous(name = "Diff",
 #                           sec.axis = dup_axis(name = element_blank()))
g <- g + theme_classic()
g <- g + theme(panel.grid.major.y = element_line(color = "black", linetype = "dotted"))+   ggplot2::labs(
    title = expression("Effect of T->G SNV on"~Delta~MFE), 
	x= expression(Delta ~MFE),
    subtitle = "Dot plot",
    caption = "source: T->G variant non 0 cosmic"
  )  

print(g)

dev.off()

#---------------------------------------------------------------

#_------------------------------------------------------------------------

temp1<-data_dfannot1%>% filter(!is.na(annot.p))%>% filter(SNP_strand=="A->G"  & which1!="no change")
%>% select(diff,name,annot.p,which1)%>% distinct()
temp1<- aggregate(diff~annot.p, temp1, FUN = median)
#temp1$diff<-temp1$diff*-1
temp1<- data.frame(temp1[order(temp1$diff), ], rank=1:nrow(temp1))

g <- ggplot(temp1, aes(y = rank, x = diff))
g <- g + geom_point(size = 3,col = "tomato3")
g <- g + scale_y_continuous(name = "Annotation", labels = temp1$annot.p, breaks = temp1$rank,
                            sec.axis = dup_axis(name = element_blank(),
                                                breaks = seq(1, nrow(temp1), (nrow(temp1)-1)/4),
                                                labels = 25 * 0:4))

#g <- g + scale_x_continuous(name = "Difference in MFE by A->G SNV",
#                            sec.axis = dup_axis(name = element_blank()))
g <- g + theme_classic()
g <- g + theme(panel.grid.major.y = element_line(color = "black", linetype = "dotted"))+   ggplot2::labs(
    title = expression("Effect of A->G SNV on"~Delta~MFE), 
	x= expression(Delta ~MFE),
    subtitle = "Dot plot",
    caption = "source: A->G variant non 0 cosmic"
  )  



print(g)

dev.off()

#---------------------------------------------------------------
temp1<-data_dfannot1%>%  filter(!is.na(annot.p))%>%
 filter(!is.na(annot.p))%>% filter(SNP_strand=="G->T" & which1!="no change")%>% select(diff,name,annot.p,which1)%>% distinct()
temp1<- aggregate(diff~annot.p, temp1, FUN = median)
#temp1$diff<-temp1$diff*-1
temp1<- data.frame(temp1[order(temp1$diff,decreasing=TRUE), ], rank=1:nrow(temp1))

g <- ggplot(temp1, aes(y = rank, x = diff))
g <- g + geom_point(size = 3,col = "darkgreen")
g <- g + scale_y_continuous(name = "Annotation", labels = temp1$annot.p, breaks = temp1$rank,
                            sec.axis = dup_axis(name = element_blank(),
                                                breaks = seq(1, nrow(temp1), (nrow(temp1)-1)/4),
                                                labels = 25 * 0:4))

#g <- g + scale_x_continuous(name = "Difference in MFE by G->T SNV",
#                            sec.axis = dup_axis(name = element_blank()))
g <- g + theme_classic()
g <- g + theme(panel.grid.major.y = element_line(color = "black", linetype = "dotted"))+   ggplot2::labs(
    title = expression("Effect of G->T SNV on"~Delta~MFE), 
	x= expression(Delta ~MFE),
    subtitle = "Dot plot",
    caption = "source: G->T variant non 0 cosmic"
  )  


print(g)

dev.off()

#---------------------------------------------------------------

#_------------------------------------------------------------------------

temp1<-data_dfannot1%>%   filter(!is.na(annot.p))%>%
filter(SNP_strand=="G->A"  & which1!="no change")%>% select(diff,name,annot.p,which1)%>% distinct()
temp1<- aggregate(diff~annot.p, temp1, FUN = median)
#temp1$diff<-temp1$diff*-1
temp1<- data.frame(temp1[order(temp1$diff,decreasing=TRUE), ], rank=1:nrow(temp1))

g <- ggplot(temp1, aes(y = rank, x = diff))
g <- g + geom_point(size = 3,col = "darkgreen")
g <- g + scale_y_continuous(name = "Annotation", labels = temp1$annot.p, breaks = temp1$rank,
                            sec.axis = dup_axis(name = element_blank(),
                                                breaks = seq(1, nrow(temp1), (nrow(temp1)-1)/4),
                                                labels = 25 * 0:4))

#g <- g + scale_x_continuous(name = "Difference in MFE by G->A SNV",
#                            sec.axis = dup_axis(name = element_blank()))
g <- g + theme_classic()
g <- g + theme(panel.grid.major.y = element_line(color = "black", linetype = "dotted"))+   ggplot2::labs(
    title = expression("Effect of G->A SNV on"~Delta~MFE), 
	x= expression(Delta ~MFE),
    subtitle = "Dot plot",
    caption = "source: G->A variant non 0 cosmic"
  )  



print(g)

dev.off()

#---------------------------------------------------------------












  ggstatsplot::theme_ggstatsplot()
#-------------------------------------------------------------------------------------------------
temp1<-data_dfannot1%>% na.omit(g4genestrand)%>%filter(diff > 10 | diff < -5)%>%
 filter((SNP_strand=="A->G") & (g4genestrand== "same" | g4genestrand=="opposite"))%>%
 select(g4genestrand,name,rel_loc,diff)%>%distinct()
temp1$g4genestrand<-as.factor(as.character(temp1$g4genestrand))
manhplot <- ggplot(temp1, aes(x = rel_loc, y = -diff, 
                                  color = g4genestrand))  + 
  geom_point(alpha = 0.75) +
  labs(title=expression("Effect of A->G SNV on"~Delta~MFE),
x = "Relative Location of SNV in a G4 region", 
       y = expression(Delta ~MFE)) + 
  theme_minimal() +
  theme( 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 60, size = 12, vjust = 0.5)
  )
manhplot
dev.off()





#------------------------------------------------
temp1<-data_dfannot1%>%filter(diff > 10 | diff < -10)%>% filter(SNP_strand=="T->G")%>% select(name,rel_loc,diff)%>%distinct()

manhplot <- ggplot(temp1, aes(x = rel_loc, y = -diff, 
                                  color = diff))  + 
  geom_point(alpha = 0.75) +
  labs(title=expression("Effect of T->G SNV on"~Delta~MFE),
x = "Relative Location of SNV in a G4 region", 
       y = expression(Delta ~MFE)) + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 60, size = 12, vjust = 0.5)
  )
manhplot
dev.off()

#-------------------------------------------------------------------------------    


data_dfannot1%>%  filter(which1=="Destabilized")%>%
 group_by(annot.p)%>%#filter(diff> quantile(diff,0.90))%>% 
summarise(n=n_distinct(name),n_genes=n_distinct(annot.gene_id),mean.num = mean(diff),
median=median(diff),q90 = quantile((diff), probs = 0.90),
mean_grt90=mean(diff[(diff)>q90]) , 
  max=min(diff),min=max(diff),q95=quantile((diff),probs=0.95),n_grt90=sum((diff)>q90))%>% distinct()%>%as.data.frame()


temp<-data_dfannot1%>%filter(which1=="Further stabilized")
foo <- pairwise.wilcox.test((temp$diff), (temp$annot.p),exact=TRUE, conf.int=TRUE, subset=1000, p.adjust.method="fdr")
foo

data_dfannot1%>%  filter(which1=="Further stabilized")%>% #select(diff,annot.gene_id,annot.p,which1)%>% distinct()%>%
 group_by(which1,annot.p)%>%
summarise(n=n_distinct(name),n_genes=n_distinct(annot.gene_id),mean.num = mean(diff),
median=median(diff),q90 = quantile(abs(diff), probs = 0.90),
mean_grt90=-mean(diff[abs(diff)>q90]) , 
  max=min(diff),min=max(diff),q95=quantile(abs(diff),probs=0.95),n_grt90=sum(abs(diff)>q90))%>% distinct()%>%as.data.frame()

 std= sqrt(var(diff)),se.num = std/ sqrt(N.num),
          
#--------

 dt<-data_dfannot1%>% filter(gquad ==  "4" | gquad== "5" |  gquad== "6")

 dt$gquad<-as.factor(dt$gquad)
 dt$which1<-as.factor(dt$which1)
dt$SNP_strand<-as.factor(dt$SNP_strand)
dt$annot.p<-as.factor(dt$annot.p)
shapiro.test(result1$residuals)

result1<-aov(diff~which1*annot.p*SNP_strand,data=dt)
library(broom)
output <- tidy(result1)
output$variationSource <- output$sumsq/sum(output$sumsq)*100
output$sig <- cut(output$p.value, breaks=c(0,  0.001,  0.01,  0.05,  0.1, 1), 
                                  labels=c( '***',  '**' , '*' , '.', ' ' ))
#---------------------
a<-data_dfannot1 %>% arrange(annot.symbol)%>%group_by(name)%>% filter(!(is.na(geneid) | annot.strand=="*"))%>% 
 summarise(n=n(),annot=paste(unique(annot.p),collapse="|"),gene.id=paste(unique(geneid),collapse="|"),genesymbol=paste(unique(symbol),collapse="|"))%>%
 arrange(desc(n))%>% as.data.frame()


data_dfannot1%>% filter(annot.p!="Intergenic" & annot.p!="lncRNA GENCODE" & !is.na(symbol))%>% head()



data_dfannot1%>% filter(annot.p!="Intergenic" & annot.p!="lncRNA GENCODE" & !is.na(symbol))%>%
filter(((V12==1 | V12==-1) & V16==0) | ((V12==0) & (V16==1 | V16==-1)))%>% summarise(gstrand=paste0(V12,":",V16),annot.p=annot.p)%>%group_by(annot.p,gstrand)%>% 
summarise(n=n())%>% dcast(annot.p~gstrand)
data_dfannot1%>% group_by(annot.strand,G_quadstrand)



c(broken_g4genes$annot.gene_id,newformed_g4genes$annot.gene_id)%>% unique()

broken_g4genes<-data_dfannot1%>%filter(annot.p=="lncRNA GENCODE")%>%filter(((V12==1) | (V12==-1)) & (V16==0))%>% filter(!is.na(annot.gene_id))%>%  select(annot.gene_id)%>% distinct()
newformed_g4genes<-data_dfannot1%>%filter(annot.p=="lncRNA GENCODE")%>% filter((V12==0)  & ((V16==-1) |(V16==1)))%>% filter(!is.na(annot.gene_id))%>%  select(annot.gene_id)%>% distinct()

 as.numeric(setdiff(as.numeric(newformed_g4genes$annot.gene_id),as.numeric(broken_g4genes$annot.gene_id)))%>% 
as.data.frame()%>% distinct()%>%  cb(.)
newformed_g4genes%>% distinct()%>% cb(.)
#select(V12,V16,name,annot.p,G_quadstrand,annot.strand,annot.symbol,typemut,g4genestrand,annot.gene_id,rel_loc,which1)%>%
 

 summarise(n=n(),annot=paste0(unique(annot.p),collapse="|"),V12,V16,gene.id=paste0(na.omit(unique(geneid)),collapse=" "),genesymbol=paste0(unique(symbol),collapse="|"),g4gene=paste0(unique(g4genestrand),collapse="|"),
which=paste0(unique(which1),collapse=""))%>% distinct()%>% arrange(desc(n))%>% as.data.frame()


data_dfannot1%>% group_by(name,G_quadstrand,geneid)%>% filter(!is.na(geneid))%>% 
 summarise(n=n(),annot=paste0(unique(annot.p),collapse="|"),gene.id=paste0(unique(geneid),collapse="|"),genesymbol=paste0(unique(symbol),collapse="|"),g4gene=paste0(unique(g4genestrand),collapse="|"),
which=paste0(unique(which1),collapse=""))%>%
 arrange(desc(n))%>% as.data.frame()%>%head()


data_dfannot1%>%
 dplyr::group_by(annot.p,g4genestrand) %>%
 dplyr::summarise(n=n()) %>%
 tidyr::spread(g4genestrand, n, fill=0) %>%
 as.data.frame()%>% cb(.)

df1<- df1%>%
   mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))

data_dfannot$V8<-gsub("\\s*\\([^\\)]+\\)","",as.character(data_dfannot$V8))
library(stringr)
data_dfannot$real_mut=gsub("^.*[||]", "", data_dfannot$V1)
data_dfannot$real_mut=gsub("&", "->", data_dfannot$real_mut)

data_dfannot$V8<-str_split(data_dfannot$V8, "\\.", simplify=T)[,1]



cb <- function(df, sep="\t", dec=",", max.size=(200*1000)){
    # Copy a data.frame to clipboard
    fwrite(df, paste0("clipboard-", formatC(max.size, format="f", digits=0)), sep=sep, row.names=FALSE, dec=dec)
  }


data_dfannot1<-left_join(data_total_2,data_dfannot[,c("V1","V5","V8","real_mut","V9","V10","V11","V12","V15","V17","V19")],by=c("name"="V1"))

transitions<-c("A->G","G->A","C->T","T->C")
transversions<-c("A->T","A->C","C->A","C->G","T->A","T->G","G->C","G->T")
data_dfannot1$typemut<-ifelse(data_dfannot1$real_mut %in% transitions,"Transition","Transversion")

temp<-split(data_dfannot1,data_dfannot1$G_quadstrand)

df_annot1 %>% 
mat<-map(temp,~xtabs(~V8+paste0(real_mut,"::"),.x))
mat<-map(temp,~xtabs(~real_mut,.x))


mat<-map(temp,~xtabs(~which1+paste0(real_mut),as.matrix(.x)))


temp<-map(temp,~setDT(.x))


df <- map(temp,~dcast(.x, SNP_strand ~ V8, length))
map(temp,~dcast(.x,V8~SNP_strand , length)[, round(.SD/Reduce(`+`, .SD)*100,2), V8])%>%
 rbindlist(id="source")%>% cb(.)


setDT(temp[[1]])[order(SNP_strand), .(SNP_strand= unique(SNP_strand),
 percentage = tabulate(SNP_strand)/.N), by = V8]

temp<-map(temp,~setDT(.x))
map(temp,~.x[ , .N, by = .(which1, real_mut) ]%>% 
dcast(which1~V8,value.var="N"))%>%rbindlist(id="source")%>% cb(.)


temp<-map(temp,setDT(.x))
map(temp,~ .x[, .N, by = .(V8,which1) ]%>% 
dcast(paste0(V8)~which1,value.var="N"))%>%rbindlist(id="source")#%>% cb(.)

map(temp,~.x[ , .N, by = .(real_mut,which1) ]%>% 
dcast(which1~real_mut,value.var="N"))%>%rbindlist(id="source")%>% cb(.)


map(temp,~.x[ , .N, by = .(V8,which1) ]%>% 
dcast(paste0(V8)~which1,value.var="N")%>% arrange(desc(Destabilized)))# cb(.)
map(temp,~.x[ , .N, by = .(left, SNP_strand) ]%>% dcast(paste0(SNP_strand,"::")~left,value.var="N"))%>% cb(.)
temp_prom<- map(temp,~.x%>% filter(V8=="promoter-TSS"))

map(temp_prom,~[ .x, .N, by = .(V8, SNP_strand,which1) ]%>% dcast(paste0(SNP_strand,"::",V8)~which1,value.var="N")#%>% cb(.)
map(temp_prom,~.x[ , .N, by = .(SNP_strand,which1) ]%>% dcast(paste0(SNP_strand,"::")~which1,value.var="N")%>% arrange(desc(Destabilized))) #cb(.)
map(temp_prom,~.x[ , .N, by = .(left, SNP_strand) ]%>% dcast(left~paste0(SNP_strand,"::"),value.var="N"))#%>% cb(.)




map(temp,~.x %>% filter(!is.na(V8))%>% filter(V8 %in% c("3' UTR","TTS","exon","promoter-TSS") &
real_mut %in% c("A->G","G->A","C->T","T->C"))%>% 
  ggplot( aes(x=rel_loc,..count..,fill=factor(SNP_strand))) +
    geom_density( alpha=0.8) +#,fill="#69b3a2", color="#e9ecef",
scale_x_continuous(breaks = seq(0,1,by=0.3))+
    ggtitle("Distribution of SNV across G4")+facet_wrap(~real_mut,scales="free")+
theme(axis.text.y=element_text(size=rel(0.8)),axis.text.x=element_text(size=rel(0.7)))) 





map(temp,~.x %>%   filter(!is.na(V8))%>%
 filter(V8 %in% c("3' UTR","TTS","exon","intron","promoter-TSS") &
SNP_strand %in% c("A->G","G->A","T->G","G->T"))%>% 
  ggplot( aes(x=rel_loc,..count..,fill=factor(typemut))) +
    geom_density( alpha=0.8) +#,fill="#69b3a2", color="#e9ecef",
scale_x_continuous(breaks = seq(0,1,by=0.3))+theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+xlab("Relative location of SNV along a G quadruplex")+
 guides(fill=guide_legend("Mutation Type"))+ 
    ggtitle("Distribution of SNV across G4")+facet_wrap(~factor(SNP_strand),scales="free")+
theme(axis.text.y=element_text(size=rel(0.8)),axis.text.x=element_text(size=rel(0.7)))) 

dev.off()


map(temp,~.x %>%   filter(!is.na(V8))%>%
 filter(V8 %in% c("3' UTR","TTS","exon","intron","promoter-TSS") &
SNP_strand %in% c("A->G","G->A","T->G","G->T"))%>% 
  ggplot( aes(x=rel_loc,..count..,fill=factor(which1))) +
    geom_density( alpha=0.8) +#,fill="#69b3a2", color="#e9ecef",
scale_x_continuous(breaks = seq(0,1,by=0.3))+theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+xlab("Relative location of SNV along a G quadruplex")+
 guides(fill=guide_legend("Mutation Type"))+ 
    ggtitle("Distribution of SNV across G4")+facet_wrap(~factor(SNP_strand),scales="free")+
theme(axis.text.y=element_text(size=rel(0.8)),axis.text.x=element_text(size=rel(0.7)))) 

dev.off()

map(temp,~.x %>%  filter(!is.na(V8)) %>%
 filter(V8 %in% c("Intergenic","exon","intron","promoter-TSS") &
SNP_strand %in% c("A->G","G->A","T->G","G->T")) %>% 
  ggplot( aes(x=rel_loc,..count..,fill=factor(which1))) +
    geom_density( alpha=0.8) +#,fill="#69b3a2", color="#e9ecef",
scale_x_continuous(breaks = seq(0,1,by=0.3))+theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+xlab("Relative location of SNV along a G quadruplex")+
 guides(fill=guide_legend("Mutation Type"))+ 
    ggtitle("Distribution of SNV across G4")+facet_wrap(V8~factor(SNP_strand),scales="free")+
theme(axis.text.y=element_text(size=rel(0.8)),axis.text.x=element_text(size=rel(0.7)))) 

dev.off()




# Make the histogram
map(temp,~.x %>% filter(!is.na(V8))%>%
  ggplot( aes(x=rel_loc,fill=which1)) +
    geom_density(aes(y=..count..), alpha=0.8) +#,fill="#69b3a2", color="#e9ecef",
scale_x_continuous(breaks = seq(0,1,by=0.1))+
    ggtitle("Distribution of SNV across G4")) #+
   # theme_ipsum()
map(temp,~.x %>% filter(!is.na(V8))%>%
  ggplot( aes(x=rel_loc,..count..,fill=V8)) +
    geom_density( alpha=0.8) +#,fill="#69b3a2", color="#e9ecef",
scale_x_continuous(breaks = seq(0,1,by=0.3))+
    ggtitle("Distribution of SNV across G4")+ facet_grid(V8~left,scales="free")+
theme(axis.text.y=element_text(size=rel(0.8)),axis.text.x=element_text(size=rel(0.7))))
 
map(temp,~.x %>% filter(!is.na(V8))%>%
  ggplot( aes(x=rel_loc,..count..,fill=V8)) +
    geom_density( alpha=0.8) +#,fill="#69b3a2", color="#e9ecef",
scale_x_continuous(breaks = seq(0,1,by=0.3))+
    ggtitle("Distribution of SNV across G4")+ facet_grid(V8~SNP_strand,scales="free")+
theme(axis.text.y=element_text(size=rel(0.8)),axis.text.x=element_text(size=rel(0.7)))) 

dev.off()
## we use this dataframe to calculate the relative position of the variant using file having start and end position of G4
##because we have used G4 based on SNVs, upstream and downstream +/- 25bp of each SNP.
# We identify G4 in each of the sequence, and based on the SNV distance from the start and end of the G4, get
# we got this using g4 parser, from names, we extract start and end position and divide by the length of the G4.
#codes used
 #./quadparser_wrapper.pl ./NCV/NCV_1_ref_g4.fasta ./NCV/NCV_1_ref_g4_quadparser.tsv

#compare the value we get next on different annotations obtained earlier.

data_dfannot[ , .N, by = .(SNP_strand, which) ]%>% dcast(SNP_strand~which,value.var="N")%>% cb(.)
data_dfannot[ , .N, by = .(V19, SNP_strand) ]%>% dcast(SNP_strand~V19,value.var="N")#%>% cb(.)
table(data_dfannot$V8,data_dfannot$SNP_strand)

temp_df<-data_dfannot1%>% filter(V19=="protein-coding")
temp_df[, sum := .N, by = V8][, prop := .N, by = list(SNP_strand, which,V8) ][, prop := prop/sum][, sum := NULL]%>% dcast(SNP_strand~which,value.var="N")
round(prop.table(table( temp_df$V8, temp_df$SNP_strand),1),2)
temp_df[ , .N, by = .(V8,real_mut) ]%>% dcast(real_mut~V8,value.var="N")%>% cb(.)

df1<-data_dfannot

##get X nucleotide context
g4_maf<-data_df
genome<-BSgenome.Hsapiens.UCSC.hg38


 g4_maf$loc_end<-g4_maf$loc_start
head(g4_maf)
g4_maf<-setDT(g4_maf)
#

# we do this because the BSgen
g4_maf<- g4_maf%>% 
       mutate(orig_refalt = case_when((G_quadstrand=="+") ~ g4_maf$SNP_strand,
                                        (G_quadstrand=="-") ~ chartr('ATGC', 'TACG',g4_maf$SNP_strand))) 


g4_maf<-separate(g4_maf,col="orig_refalt",into=c("ref","alt"),remove=FALSE)


#g4_maf$chr<-paste0("chr",g4_maf$chr)

  
g4_maf$loc_end<-g4_maf$loc_start
#g4_maf$chr<-paste0("chr",g4_maf$chr)
g4_maf<-g4_maf[((nchar(g4_maf$ref)==1) | (nchar(g4_maf$alt)==1)),]
seq_g4_1 <- BSgenome::getSeq(genome, names = g4_maf$left,
                             start = (g4_maf$loc_start-1),
                             end = (g4_maf$loc_start-1))

seq_g4_2 <- getSeq(genome, names = g4_maf$left,
                   start = (g4_maf$loc_start+1),
                   end =  (g4_maf$loc_start+1))

seq_g4_full <- getSeq(genome, names = g4_maf$left,
                      start = (g4_maf$loc_start-1),
                      end = g4_maf$loc_start+1)


head(as.data.frame(seq_g4_1))
head(as.data.frame(seq_g4_2))
head(as.data.frame(seq_g4_full))
g4_seq<-(seq_g4_full)
g4_maf <- g4_maf %>% mutate(G4_id = row_number())
#g4_maf_s<-g4_maf[,c(1,2,3,9,14,15,16,17,18,19,21,22,43,44,45,126,127,128,129,130)]
#g4_maf_s<-g4_maf[,c(1,2,3,7,8,12,13,14,15,16,17,18,19,20,22,23,24,25,26,27,121,122,123,124,125)]
g4_maf_s<-g4_maf#[,c(1,2,3,4,5,6,9,10,11,12,15,16,18,19,20)]
seq_left<-as.data.frame(seq_g4_1)
colnames(seq_left)<-"leftseq"
#last letter of the left seq is the position in question
seq_right<-as.data.frame(seq_g4_2)


revString <- function(text){
  paste(rev(unlist(strsplit(text,NULL))),collapse="")
}
 

rm(seq_g4_1,seq_g4_2)
colnames(seq_right)<-"rightseq"
seq_g4_full<-as.data.frame(seq_g4_full)
colnames(seq_g4_full)<-"seq_g4_full"
seq<-cbind(as.data.frame(seq_left),as.data.frame(g4_maf_s),as.data.frame(seq_right),as.data.frame(seq_g4_full))

seq$G4_ref<-paste0(seq$leftseq,seq$ref,seq$rightseq)
seq$G4_alt<-paste0(seq$leftseq,seq$alt,seq$rightseq)
                                                                             
all_1<-seq


all_2<- all_1 %>% #rowwise()%>%
       mutate(G4_ref = case_when((G_quadstrand=="+") ~ all_1$G4_ref,
                                        (G_quadstrand=="-") ~ intToUtf8(rev(utf8ToInt(G4_ref)))),
		G4_alt=case_when((G_quadstrand=="+") ~ all_1$G4_alt,
                                        (G_quadstrand=="-") ~ intToUtf8(rev(utf8ToInt(G4_alt)))))


all_2<-all_1%>% rowwise()%>% 
mutate(G4_ref=ifelse(all_1$G_quadstrand=="-",intToUtf8(rev(utf8ToInt(G4_ref))),
G4_ref),
G4_alt=ifelse(all_1$G_quadrstrand=="-",intToUtf8(rev(utf8ToInt(G4_alt))),G4_alt))

reversed_string <- intToUtf8(rev(utf8ToInt(G4_ref)))

#tricosmic<-all_2
pencosmic<-all_2
save(tricosmic,pencosmic,file="/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/flankSNVcontent_G4cosmic.RData")

#load("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/flankSNVcontent_G4cosmic.RData")
tricosmic<-setDT(tricosmic)
aa<-tricosmic%>% mutate(aref=str_count(G4_ref,pattern="G|C"),aalt=str_count(G4_alt,pattern="G|C"))%>%
filter(aref>1 | aalt>1)#%>%# filter(G4_ref~ !grepl('G|C') | G4_alt ~ !grepl('G|C') )
aa[ , .N, by = .(G4_ref, G4_alt) ]%>%# mutate(aref=str_count(G4_ref,pattern="G|C"),aalt=str_count(G4_alt,pattern="G|C"))%>%
2#filter(aref>1 | aalt>1)%>%
dcast(G4_ref~G4_alt,value.var="N")#%>% cb(.)




data_df<-setkey(data_df, left, start, end)
annotation<-setkey(annotation, V1, V2, V3)

resultTable<-foverlaps(data_df,annotation,nomatch=0)
resultTable <- map(data,~foverlaps(.x, annotation, nomatch = 0))

df.g <- group_split(df.h, SNP)
map_dfr(df.g, function(x) tidy(lm(value ~ wind + temp, data=x)))
g4_refalt_rnafold%>%group_by(loc,SNP,which,region)%>% summarise(p=paste(file[order(file)],collapse="|"),n=n_distinct(file))%>%
  as.data.frame()%>% arrange(n)%>% nrow()


#g4_refalt_rnafold<-  merge(RNAFOLD_NCV,RNAFOLD_coding,by=c("name","MFE.x","MFEadj.x","EFE.x","MFE.y", "MFEadj.y","EFE.y" ))
g4_refalt_rnafold<-rbind(RNAFOLD_NCV,RNAFOLD_coding)%>%as.data.frame()%>%arrange(file)%>% 
  group_by(name,MFE.x,MFEadj.x, CDE.x ,MEAFE.x, MEA.x, ED.x, MFEadj.GC.x, MEAFEadj.x,
           EDadj.x,  CFE.x ,EFE.x,MFE.y, MFEadj.y,EFE.y, CDE.y, MEAFE.y, MEA.y ,ED.y,
           MFEadj.GC.y, MEAFEadj.y ,   EDadj.y , CFE.y )%>%
  summarise(file=paste0(unique(file),collapse = "|"))%>% distinct()%>% as.data.frame()



ref_RNAFOLD_coding<-setDT(ref_RNAFOLD_coding)
coding2<-setDT(coding2)
coding_merged<-left_join(ref_RNAFOLD_coding,coding2,by="name")

coding_merged<-coding_merged%>% select(-c(G4_ref,G4_alt,start,end,left))



#for clinvar:
get total for most variables:
#annotate the CLINVAR AND COSMIC whole database:
#awk 'BEGIN{FS=OFS="\t"} {print "chr"$2,$3,$4,NR==1 ? "ID" : "\""NR-1"\""$1,$5,$7}' COSMIC_merged1.bed 
#    awk 'BEGIN{FS=OFS="\t"} {$2="chr"$2; print $2,$3,$4,"identifier"$1,$5,$7}' | head

#awk 'OFS="\t" { $1="chr"$1; print}' in.gtf


COSMIC_all<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/COSMIC_all_SNV_annotated.bed")


COSMIC_all$Annotation<-gsub("\\s*\\([^\\)]+\\)","",as.character(COSMIC_all$Annotation))
library(stringr)

COSMIC_all$Annotation<-str_split(COSMIC_all$Annotation, "\\.", simplify=T)[,1]

 tempcosmic<-fread("COSMIC_merged1.bed")
colnames(COSMIC_all)[1]<-"ID"
COSMIC_all<-left_join(COSMIC_all,tempcosmic[,c("V4","V6","V7")],by=c("ID"="V4"))
COSMIC_all<-COSMIC_all%>% filter(nchar(V6)==1 & nchar(V7)==1)

total_count<-COSMIC_all[ , .N, by = .(V6,V7,Annotation) ]%>% 
dcast(paste0(V6,"->",V7)~Annotation,value.var="N")#%>% cb(.)

baseR.sbst.rssgn   <- function(x) { x[is.na(x)] <- 0; x }
total_count<-baseR.sbst.rssgn(total_count)
cbind(id = total_count[, 1], round(total_count[, -1]*100/rowSums(total_count[, -1]),2))%>% cb()

temp_count<-map(temp,~.x[ , .N, by = .(SNP_strand,V8) ]%>% dcast(SNP_strand~V8,value.var="N"))#%>% 
#rbindlist(source="id")%>% cb(.)
#temp_count<-map(temp_count,~.x[,-c(2)])
temp_count<-map(temp_count,~baseR.sbst.rssgn(.x))

map(temp_count,~cbind(id = .x[, 1], round(.x[, -1]*100/rowSums(.x[, -1]),2)))%>% rbindlist(idcol="id")%>%cb(.)
temp_count%>% rbindlist(idcol="id")%>% cb()







##trying to calculate SNVs in the flanking region of G4

bedtools closest  -d -t all  -a variants_aroundG4.cosmic.bg -b G4_aroundSNVcosmic_all.sortedmerged.bed > variants_aroundG4.cosmic.withG4closest.bed

bedtools intersect  -loj -a variants_aroundG4.cosmic.bg -b G4_aroundSNVcosmic_all.sortedmerged.bed > variants_aroundG4.cosmic.withG4.bed


flankG4SNV<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/variants_aroundG4.cosmic.withG4closest.bed")

flankG4SNV0<-flankG4SNV%>% filter(V10==0)
flankG4SNVflanknot0<-flankG4SNV%>% filter(V10!=0)
flankG4SNVflanknot0$distgroup<-cut(flankG4SNVflanknot0$V10,seq(from=-2000,to=2000,by=100),right=FALSE,include.lowest=FALSE)
flankG4SNV0$distgroup<-"0"
flankSNV<-rbind(flankG4SNV0,flankG4SNVflanknot0)
rm(flankG4SNV0,flankG4SNVflanknot0)
rm(flankG4SNV)
flankSNV$distgroup<-factor(flankSNV$distgroup)

flankSNV$distgroup1<-ordered(flankSNV$distgroup,levels=c("[-2e+03,-1.9e+03)",
   "[-1.9e+03,-1.8e+03)", "[-1.8e+03,-1.7e+03)",
"[-1.7e+03,-1.6e+03)" ,"[-1.6e+03,-1.5e+03)", "[-1.5e+03,-1.4e+03)",
 "[-1.4e+03,-1.3e+03)", "[-1.3e+03,-1.2e+03)", "[-1.2e+03,-1.1e+03)",
 "[-1.1e+03,-1e+03)",   "[-1e+03,-900)",       "[-900,-800)",
"[-800,-700)"   , "[-700,-600)",  "[-600,-500)",
"[-500,-400)" , "[-400,-300)",    "[-300,-200)",
 "[-200,-100)" , "[-100,0)"   ,  "0" , "[0,100)",
 "[100,200)", "[200,300)" , "[300,400)",
 "[400,500)",  "[500,600)" , "[600,700)",
"[700,800)" , "[800,900)"  , "[900,1e+03)",
 "[1e+03,1.1e+03)" , "[1.1e+03,1.2e+03)"   ,"[1.2e+03,1.3e+03)",
 "[1.3e+03,1.4e+03)", "[1.4e+03,1.5e+03)"   ,"[1.5e+03,1.6e+03)",
 "[1.6e+03,1.7e+03)" ,  "[1.7e+03,1.8e+03)"  , "[1.8e+03,1.9e+03)",
 "[1.9e+03,2e+03]" ,"[1.9e+03,2e+03)"))

temp<-split(flankSNV,flankSNV$V9)

remove 
#    summarise(across(where(is.numeric), 
#          ~ mean(.[!. %in% range(.)])), .groups = 'drop')


map(temp,~ggplot(.x, aes(x = distgroup, y = log(V4),fill=distgroup)) +   geom_boxplot() +
    xlab("class") +
    theme(legend.position="none") +
    xlab("SNV in/flanking(+/-2000bp) G4"))



tgc<-map(temp,~.x %>%
  dplyr::group_by(distgroup1) %>%
  summarise(mean.count = mean(V4, na.rm = TRUE),
            sd.count = sd(V4, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se.count = sd.count / sqrt(n.count),
         lower.ci.count = mean.count - qt(1 - (0.05 / 2), n.count - 1) * se.count ,
         upper.ci.count = mean.count + qt(1 - (0.05 / 2), n.count - 1) * se.count )%>% as.data.frame())#%>% 
#rbindlist(idcol="source")%>% as.data.frame()%>%fwrite(.,file="G4_comsmic_Flank_final.tsv",sep="\t")

pd <- position_dodge(0.1) # move them .05 to the left and right
pdf(file= "/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/pg4_COSMIC_flank_count.pdf",width=10,height=7)
names(tgc)<-c("positive strand","Negative Strand")
map(tgc,~ggplot(.x, aes(x=distgroup1, y=mean.count)) + 
    geom_errorbar(aes(ymin=mean.count-lower.ci.count, ymax=mean.count+upper.ci.count), colour="black", width=.2, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3)+scale_colour_Publication()+
theme_bw()+ xlab("Distance from G quadruplex divided into equal 100 bp window") + ylab("Counts of SNV")+
# theme_Publication()+
theme(axis.text.x = element_text(angle = 30, hjust = 1), 
legend.title=element_blank(), legend.position="bottom", legend.key = element_rect(color = 'white')))


dev.off()

df<-flankSNV%>% group_by(distgroup)%>% count()
df<-flankSNV%>% group_by(distgroup)%>% summarise(mean1=mean(V4),sum=sum(V4),median=median(V4),sd=sd(V4))

summary(flankSNV)
one.way <- aov(V4 ~ distgroup1, data = flankSNV)
interaction <- aov(V4 ~ distgroup1*V9, data = flankSNV)

df$distgroup<-ordered(df$distgroup,levels=c("[-2e+03,-1.9e+03)",
   "[-1.9e+03,-1.8e+03)", "[-1.8e+03,-1.7e+03)",
"[-1.7e+03,-1.6e+03)" ,"[-1.6e+03,-1.5e+03)", "[-1.5e+03,-1.4e+03)",
 "[-1.4e+03,-1.3e+03)", "[-1.3e+03,-1.2e+03)", "[-1.2e+03,-1.1e+03)",
 "[-1.1e+03,-1e+03)",   "[-1e+03,-900)",       "[-900,-800)",
"[-800,-700)"   , "[-700,-600)",  "[-600,-500)",
"[-500,-400)" , "[-400,-300)",    "[-300,-200)",
 "[-200,-100)" , "[-100,0)"   ,  "0" , "[0,100)",
 "[100,200)", "[200,300)" , "[300,400)",
 "[400,500)",  "[500,600)" , "[600,700)",
"[700,800)" , "[800,900)"  , "[900,1e+03)",
 "[1e+03,1.1e+03)" , "[1.1e+03,1.2e+03)"   ,"[1.2e+03,1.3e+03)",
 "[1.3e+03,1.4e+03)", "[1.4e+03,1.5e+03)"   ,"[1.5e+03,1.6e+03)",
 "[1.6e+03,1.7e+03)" ,  "[1.7e+03,1.8e+03)"  , "[1.8e+03,1.9e+03)",
 "[1.9e+03,2e+03]" ,"[1.9e+03,2e+03)"))


ggplot(df, aes(x=1,y=distgroup,weights=mean1)) + geom_boxplot()

flankG4SNV_combined<-left_join(flankG4SNV,
data_dfannot1[,c("name","loc_start","affectg4","numberofg4","V8",
"rel_loc","which1","V10","V15","V19")],by=c("V8"="name"))


flankG4SNV_combined %>% filter(affectg4=="affect_g4")%>% 
group_by(V8,V10.x)%>%mutate(V4=as.numeric(V4))%>%
 summarise(mean1=mean(V4),median1=median(V4),sd1=sd(V4))%>% group_by(V10.x)%>%
  summarise(mean1=mean(mean1),median1=median(median1),sd1=sd(sd1))%>%distinct()%>%as.data.frame()%>% head()


flankG4SNV%>% (dist=V5-V2)%>% group_by(V8)


tricosmic
 temp<-tricosmic%>% group_by(real_mut)%>% 
mutate(countmut=n())%>% group_by(G4_ref)%>% 
mutate(refcount=n())%>% group_by(G4_alt)%>% 
 summarise(real_mut,G_quadstrand,G4_ref,refcount,countmut,altcount=n())%>%
distinct()%>% as.data.frame()
















#!/usr/bin/bash -l 
cat $@ | \ 
    sed 's/"//g' | `# remove quotes` \
    awk 'function notempty_or(x,y) {return x?x:y} \
         BEGIN{FS="\t"} \
         { \
            if($3!~"gene") { next }; \
            split($9, a, "; "); \
            for (i in a) { \
                split(a[i], j, " "); \
                attr[ j[1] ] = j[2]; \
            };\
            gene_id = notempty_or(attr["gene_id"],"NA"); \
            gene_name = notempty_or(attr["gene_name"],gene_id); \
            print gene_id "\t" gene_name ; \
            attr["gene_id"]=""; \
            attr["gene_name"]=""; \
         }' | \ 
    sed 's/;$//' |  `# remove terminal semicolons` \
    sed --regexp-extended 's/;(\s)/\1/' | `# remove semicolon before space character, but keep the space` \
    sort