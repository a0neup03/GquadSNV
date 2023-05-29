#------------------------------------------------------------------------------------------------------------------------------------------------------

library(annotatr)
library(dplyr)
library(ggplot2)
library(purrr)
library(ggtext)
library(data.table)
#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#setwd("C:/Users/a0neup03/OneDrive - University of Louisville/Documents/SNV_paper/g4snv/")
#load("C:/Users/a0neup03/OneDrive - University of Louisville/Documents/SNV_paper/g4snv/full_clinvar_G4_duplicatesremoved.Rdata")
a<-load("C:/Users/a0neup03/OneDrive - University of Louisville/Documents/SNV_paper/g4snv/Clinvar_G4.Rdata")

head(data_total_2)

#####
#------------------------------------------------------------------------------------------------------------------------------------------------------

#this is the ultimate clinvar for G4 variants without gene annotations. here the rows are repeated. 
#so to count using data_total_2 use any specific columns and non repeated rows.
#starting with how many variants affect the G4 directly or are present just upstream or downstream of the G4
x<-data_total_2%>% distinct()
table(x$affectg4,x$G_quadstrand)
x$affectg4<-factor(x$affectg4,levels=c("affect_g4","mut_upstream","mut_downstream"))
x <- x[order(x$affectg4),]
x1<-x[!duplicated(x[,c('name', 'G_quadstrand')]),]
table(x1$affectg4,x1$G_quadstrand)

#from this we select the g4 affecting directly the variants.#this is what we use now.
x2<-x1%>% filter(affectg4=="affect_g4")

# Count rows with NA in ref and non NA in ALT based on Quadparser
sum(rowSums(is.na(x2[c("g4seq.x")])) == 1 & !is.na(x2$g4seq.y))    #401
# Count rows with non-NA in alt and NA in Ref
sum(rowSums(is.na(x2[c( "g4seq.y")])) == 1 & !is.na(x2$g4seq.x))   #1186
#with G4 unchanged by mutation
sum(rowSums(!is.na(x2[c( "g4seq.y")])) == 1 & !is.na(x2$g4seq.x)) #3412
#check if values are correct
table(lapply(x2[c("g4seq.x","g4seq.y")],is.na))

sum(rowSums(is.na(x2[c( "g4seq.y")])) == 1 & !is.na(x2$g4seq.x) & x2$V10==1)   #577
sum(rowSums(!is.na(x2[c( "g4seq.y")])) == 1 & !is.na(x2$g4seq.x) & x2$V10==1 & x2$V16==1)  #1697
sum(rowSums(is.na(x2[c( "g4seq.x")])) == 1 & !is.na(x2$g4seq.y) & x2$V16==1)   #200

#using G4hunter scores
table(x2$V10,x2$V16)


x2%>% group_by(name,G4_id,left,right,G4_ref,G4_alt,V9,V10,V11,V12,V13,V14,V15,V16,G4_refrc,G4_altrc,start,end,REF,ALT,region,loc,SNP,MFE.x,MFE.y,CDE.x,CDE.y,CFE.x,CFE.y,SNP_strand,G_quadstrand)%>%
  summarise(n=n())%>% arrange(desc(n))%>%nrow()# head()%>%as.data.frame()
#4999 variants total in CLINVAR

#======================================================================================================================================================
#these are all the variants identified. Now based on these, we intersect the known experimental locations using GSM3003539_w15.hits.K.bed
#======================================================================================================================================================
#we get the hg19 files for G4 sequences using GSM3003539 and convert to hg38
chain_file=  "C:/Users/a0neup03/OneDrive - University of Louisville/Documents/hg38/hg19ToHg38.over.chain.gz"

chr <- c(1:22,"X","Y","MT")

library(rtracklayer)
library(GenomicRanges)
chain_file=  "C:/Users/a0neup03/OneDrive - University of Louisville/Documents/hg38/hg19ToHg38.over.chain"
chainObject <-import.chain(chain_file)
GSM<-makeGRangesFromDataFrame(fread("C:/Users/a0neup03/Documents/G4files/GSM3003539_w15.hits.K.bed",col.names = c("chr","start","end")))

GSM_hg38 <- as.data.frame(liftOver(GSM, chainObject))

#------------------------------------------------------------------------------------------------------------------------------------------------------
#======================================================================================================================================================

bedA<-x2%>% filter(rel_loc>=0& rel_loc<=1)%>%dplyr::select(left,start,end,name,G_quadstrand)%>% distinct()%>% 
  reframe(CHROM=left,START=start,STOP=end,strand=G_quadstrand,id=name)%>% distinct()%>%
  makeGRangesFromDataFrame(.,keep.extra.columns = TRUE)

GSM_hg38<-GSM_hg38%>% dplyr::select(group,seqnames,start,end,width,strand)%>%
  reframe(group=group,CHROM=seqnames,START=start,STOP=end,strand=strand,width=width)%>%  distinct()

GSM_hg38<-makeGRangesFromDataFrame(GSM_hg38,keep.extra.columns = TRUE)

# subset f1 features that overlap with f2
result <- bedA[countOverlaps(bedA, GSM_hg38, maxgap=-1) > 0]
G4_var_GSM<-as.data.frame(result)


#save(G4_var_GSM,file="CLINVAR_geo_expr_g4.RData")

G4_expr<-left_join(G4_var_GSM,x2,by=c("id"="name"))



#======================================================================================================================================================
#======================================================================================================================================================
#======================================================================================================================================================

library(dplyr)
require(Biostrings) 
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
length(Hsapiens)

#seqnames(Hsapiens)
genome<-BSgenome.Hsapiens.UCSC.hg38

#g4_maf<-as.data.frame(g4_var)
#g4_maf<-as.data.frame(g4_var)
#colnames(g4_maf)<-c("name","chr","start","end","G4score","chr1","loc_start","var_id1","ref","alt","gene","strand","var_id","sth","transcript","SNP","sequence")
#colnames(G4_expr)<-c("name","chr","start","end","loc_end","ref","alt","SNV","MFE_ref","MFE_alt","which.ED","which","loc_start",
 #                   "diff_start","diff_end","length","ldiff_start","ldiff_end")

#g4_maf<-g4_maf%>% group_by(name,loc_start,chr,start,end,ref,alt,SNV)  %>% summarise(n=n())%>% arrange(desc(n))

#check if bilallelic or multinucleotide variant present 
#G4_expr<-G4_expr[((nchar(G4_expr$ref)==1) | (nchar(G4_expr$alt)==1)),]
G4_expr$right<-as.numeric(G4_expr$right)
G4_expr$chr<-as.character(G4_expr$seqnames)

seq_g4_1 <- BSgenome::getSeq(genome, names = G4_expr$seqnames,
                             start = as.numeric(G4_expr$right-1),
                             end = as.numeric(G4_expr$right-1))

seq_g4_2 <-  BSgenome::getSeq(genome, names = G4_expr$seqnames,
                   start = as.numeric(G4_expr$right+1),
                   end =  as.numeric(G4_expr$right+1))

seq_g4_full <- getSeq(genome, names = G4_expr$seqnames,
                      start = as.numeric(G4_expr$right-2),
                      end = as.numeric(G4_expr$right+2))


head(as.data.frame(seq_g4_1))
head(as.data.frame(seq_g4_2))
head(as.data.frame(seq_g4_full))
g4_seq<-(seq_g4_full)
g4_maf <- G4_expr %>% mutate(G4_id = row_number())
#here we reverse REF and ALT for negative strands

#g4_maf_s<-g4_maf[,c(1,2,3,9,14,15,16,17,18,19,21,22,43,44,45,126,127,128,129,130)]
#g4_maf_s<-g4_maf[,c(1,2,3,7,8,12,13,14,15,16,17,18,19,20,22,23,24,25,26,27,121,122,123,124,125)]
g4_maf_s<-g4_maf#[,c(1,2,3,4,5,6,9,10,11,12,15,16,18,19,20)]
seq_left<-as.data.frame(seq_g4_1)
colnames(seq_left)<-"leftseq"

#last letter of the left seq is the position in question
seq_right<-as.data.frame(seq_g4_2)

rm(seq_g4_1,seq_g4_2)
colnames(seq_right)<-"rightseq"
seq_g4_full<-as.data.frame(seq_g4_full)
colnames(seq_g4_full)<-"seq_g4_full"

seq<-cbind(as.data.frame(seq_left),as.data.frame(g4_maf_s),as.data.frame(seq_right),as.data.frame(seq_g4_full))

seq$G4_ref<-paste0(seq$leftseq,seq$REF,seq$rightseq)
seq$G4_alt<-paste0(seq$leftseq,seq$ALT,seq$rightseq)

library(stringi)
#G4_refrc abd G4_altrc are both 5' to 3'. Even if G quadruplex is in the negative strand, based on the mutations we take G4 from 5' to 3'
all_1<-seq%>% mutate(mref=ifelse(strand=="+",substr(G4_refrc, 29, 31),substr(G4_refrc, 31, 33)),
                                 malt=ifelse(strand=="+",substr(G4_altrc, 29, 31),substr(G4_altrc, 31, 33)))%>%
  mutate(G4_ref=ifelse(strand=="-",stri_reverse(chartr('ATGC', 'TACG', G4_ref)),G4_ref),
         G4_alt=ifelse(strand=="-",stri_reverse(chartr('ATGC', 'TACG', G4_alt)),G4_alt))#
                                 
  
#we can check if the sequences we are comparing is correct or not
all_1%>%
  filter(strand=="+")%>% select(id,strand,mref,malt,seq_g4_full,G4_ref,G4_alt,G4_refrc,G4_altrc)%>% head()

all_1%>%
  filter(strand=="-")%>% select(id,strand,mref,malt,seq_g4_full,G4_ref,G4_alt,G4_refrc,G4_altrc)%>% rowwise()%>%
  head()%>% as.data.frame()
 

m<-table(paste0(all_1$mref,":",all_1$malt),all_1$G_quadstrand)%>% as.data.frame()%>% arrange(desc(Freq))#%>% write.table(., "clipboard-16384", sep="\t", row.names=FALSE,quote=FALSE)

m%>% head()

m<-table(as.character(paste0(all_1$G4_ref,"->",all_1$G4_alt)),all_1$SNP_strand,all_1$G_quadstrand)%>% as.data.frame()%>% arrange(desc(Freq))#%>% spread(Var3,Freq)%>% 
#  arrange(desc(`+`))%>% as.data.frame()%>%head(20)

m<-table(as.character(paste0(all_1$G4_ref,"->",all_1$G4_alt)),all_1$SNP_strand,all_1$chr)%>% as.data.frame()%>% arrange(desc(Freq))#%>% spread(Var3,Freq)%>% 

m %>% 
  group_by(Var3) %>% 
  distinct(Var1, Var2, Var3) %>%
  ungroup %>%
  summarise(pval = chisq.test(Var2, Var1)$p.value)


table(paste0((all_1$G4_ref),"->", (all_1$G4_alt)),all_1$SNP_strand,all_1$G_quadstrand)%>% as.data.frame()%>% arrange(desc(Freq))%>% spread(Var3,Freq)%>% 
  arrange(desc(`-`))%>% as.data.frame()%>%head(20)

table(as.character(paste0(all_1$mref,"->",all_1$malt)),all_1$SNP_strand,all_1$G_quadstrand)%>% as.data.frame()%>% arrange(desc(Freq))%>% spread(Var3,Freq)%>% 
  arrange(desc(`+`))%>% head(20)
x<-all_1%>% group_by(mref,malt,SNP_strand,G_quadstrand)%>% summarise(Freq=n())%>% 
  group_by(SNP_strand)%>% mutate(SNP_n=sum(Freq))%>% spread(G_quadstrand,Freq,fill = 0)%>% mutate(neg_s=(`-`*100)/SNP_n,pos_s=(`+`*100)/SNP_n)%>%
arrange(desc(`-`))%>% distinct()#%>%filter(SNP_n>50)%>%head(20)
  

#plot count of variants by strands presence in G4
all_1%>% group_by(mref,malt,SNP_strand,G_quadstrand)%>% summarise(Freq=n())%>%
group_by(SNP_strand)%>% mutate(SNP_n=sum(Freq),gg=paste0(mref,"->",malt))%>%mutate(s=(Freq*100)/SNP_n)%>%
         ggplot(.,aes(x = SNP_strand,  y = Freq, fill = SNP_strand))+# Grouped barplot using ggplot2
 xlab("SNV")+ ylab("Frequency")+
geom_bar(stat = "identity",
               position = "dodge")+facet_wrap(~G_quadstrand)+  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#plot count of variants by presence/absence of G4
all_1%>%  
  group_by(g4hunter)%>% mutate(SNP_n=sum(Freq),gg=paste0(mref,"->",malt))%>%mutate(s=(Freq*100)/SNP_n)%>%
  ggplot(.,                                      # Grouped barplot using ggplot2
         aes(x = SNP_strand,
             y = Freq,
             fill = SNP_strand)) + xlab("SNV")+ ylab("Frequency")+
  geom_bar(stat = "identity",
           position = "dodge")+facet_wrap(~g4hunter)+  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
#→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→

all_1%>% filter(!(V12==-1 & V16==1))%>% mutate(g4hunter=case_when((V12==-1 | V12==1) & (V16==1 | V16==-1)~"G4 Present",
                                               ((V12==-1 | V12==1) & (V16==0))~"G4 loss", 
                                               ((V12==0) & (V16==1 | V16==-1))~"G4 Gain",
                                                         TRUE~ "sth"))%>%
  group_by(mref,malt,SNP_strand,g4hunter)%>% summarise(Freq=n())%>% 
  group_by(g4hunter)%>% mutate(SNP_n=sum(Freq),gg=paste0(mref,"→",malt))%>%mutate(s=(Freq*100)/SNP_n)%>%
ggplot(.,                                      # Grouped barplot using ggplot2
       aes(x = SNP_strand,
           y = Freq,
           fill = SNP_strand)) + xlab("SNV")+ ylab("Frequency")+
  geom_bar(stat = "identity",
           position = "dodge")+facet_wrap(~g4hunter)+  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

#→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→
transitions<-c("A->G","G->A","C->T","T->C")
transversions<-c("A->T","A->C","C->A","C->G","T->A","T->G","G->C","G->T")

all_1$typemut<-ifelse(all_1$SNP %in% transitions,"Transition","Transversion")


all_1%>% group_by(typemut,G_quadstrand)%>% summarise(Freq=n())%>% ungroup()%>% distinct()%>%
  group_by(typemut)%>% mutate(SNP_n=sum(Freq))%>%mutate(s=(Freq*100)/SNP_n)%>%
  ggplot(.,aes(x = typemut,  y = Freq, fill = typemut))+# Grouped barplot using ggplot2
  xlab("SNV")+ ylab("Frequency")+
  geom_bar(stat = "identity",
           position = "dodge")+facet_wrap(~G_quadstrand)+  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


#→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→

all_1%>% group_by(mref,malt,typemut,G_quadstrand)%>% summarise(Freq=n())%>%
  group_by(typemut)%>% mutate(SNP_n=sum(Freq),gg=paste0(mref,"→",malt))%>%mutate(s=(Freq*100)/SNP_n)%>%filter(Freq>10)%>%
  ggplot(.,aes(x = reorder(as.character(gg),-Freq),  y = Freq, fill = typemut))+# Grouped barplot using ggplot2
  xlab("SNV")+ ylab("Frequency")+
  geom_bar(stat = "identity",
           position = "dodge")+facet_wrap(~typemut,scales="free")+  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+theme(axis.text=element_text(size=12),
                                                                               axis.title=element_text(size=14,face="bold"))


#dont waste time on this
#this becomes if we want to check for negative stranded G4 occuring 3' direction
table(paste0(stri_reverse(all_1$mref),"→", stri_reverse(all_1$malt)),all_1$SNP_strand,all_1$G_quadstrand)%>% as.data.frame()%>% arrange(desc(Freq))%>% spread(Var3,Freq)%>% 
  arrange(desc(`-`))%>% as.data.frame()%>%head(20)
#-------------------------------------------------------------------------------------
#some basic counts of G4 in clinvar based on type of mutation, effect of strand, strand etc.

tabletemp<-
  all_1%>% dplyr::select(SNP_strand,V9,V12,V16,V13,V16,id,typemut,which,rel_loc)%>% distinct()%>%
  separate(SNP_strand,into=c("ref","alt"),sep="->")

as.data.frame.matrix(table(tabletemp$ref,tabletemp$alt))%>% 
  fwrite(.,file="CLINVAR_SNV_strand.tsv",sep="\t")
as.data.frame.matrix(round(prop.table(table(tabletemp$ref,tabletemp$alt))*100,2))%>% 
  fwrite(.,file="CLINVAR_SNV_strand_prop.tsv",sep="\t",row.names=TRUE)
table(tabletemp$typemut,tabletemp$which)

as.data.frame.matrix(table(tabletemp$ref,tabletemp$alt))%>% 
  fwrite(.,file="CLINVAR_SNV_type_affect_count.tsv",sep="\t")

as.data.frame.matrix(round(prop.table(table(tabletemp$typemut,tabletemp$which))*100,2))%>% 
  fwrite(.,file="CLINVAR_SNV_type_affect_perc.tsv",sep="\t",row.names = TRUE) 
as.data.frame.matrix(table(tabletemp$typemut,tabletemp$which))%>% 
  fwrite(.,file="CLINVAR_SNV_type_affect_count.tsv",sep="\t",row.names = TRUE) 


as.data.frame.matrix(table(paste0(tabletemp$V12),tabletemp$V16))

round(tabletemp %>% group_by(V12,V16)%>% 
        summarise(mn=n(),V9ref=mean(V9),V9sd=sd(V9),
                  V13alt=mean(V13),sdalt=sd(V13))%>% as.data.frame(),3)%>%
fwrite(.,file="CLINVAR_g4hunter_count.tsv",sep="\t") 


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


m<-all_1%>%filter(!(is.na(G4_ref)| is.na(G4_alt)))%>%              
  dplyr::mutate(V1factor=case_when( grepl('^[ACTG][ATC][ACTG]$',G4_ref) | 
                                      grepl('^[ACGT][ACGT][ACT]$',G4_ref) |
                                      grepl('^[ACT][ACGT][AGCT]$',G4_ref) ~"loop",
                                    grepl('^[G][G][G]$',G4_ref) ~"G4")) %>%rowwise()%>% 
  mutate(V2factor = case_when(  grepl('^G[ATC]G$',G4_alt) | 
                                  grepl('^[ACT][ACGT][ACGT]$',G4_alt) |
                                  grepl('^[ACGT][ACGT][ACT]$',G4_alt) ~"loop",
                                grepl('^[G][G][G]$',G4_alt) ~"G4")) %>% as.data.frame()

#all_1<-setDT(all_1)
#all_1[, V1factor := ifelse(str_detect(G4_ref, "GGG"), "G4", "loop")]
#all_1[, V2factor := ifelse(str_detect(G4_alt, "GGG"), "G4", "loop")]


as.data.frame.matrix(round(prop.table(table(m$typemut,paste0(m$V1factor,":",m$V2factor)),1)*100,2))%>%
  fwrite(.,file="CLINVAR_typemut3mer_prop.tsv",sep="\t",row.names=T) 
as.data.frame.matrix((table(m$typemut,paste0(m$V1factor,":",m$V2factor))))%>%
  fwrite(.,file="CLINVAR_typemut3mer_count.tsv",sep="\t",row.names=T) 

as.data.frame.matrix((table(m$which,paste0(m$V1factor,":",m$V2factor))))%>%
  fwrite(.,file="CLINVAR_which3mer_count.tsv",sep="\t",row.names = T) 
as.data.frame.matrix(round(prop.table(table(m$which,paste0(m$V1factor,":",m$V2factor)),1)*100,2))%>%
  fwrite(.,file="CLINVAR_which3mer_prop.tsv",sep="\t",row.names=T) 

#-----------------------------------------------------------------------------------------------------------------------
#for this requires annotation information present in clinvar_annot
#we use clinvar_annot.l generated below using annotatr and separate annotation for DNA and RNA
#clinvar_annot.l

tabletemp<-purrr::map(clinvar_annot.l,~.x%>%
  dplyr::select(SNP_strand,g4genestrand,annot.p,id,typemut,which,rel_loc,V1factor,V2factor)%>% distinct())#%>%
#separate(SNP_strand,into=c("ref","alt"),sep="->"))

purrr::map(tabletemp,~round(prop.table(table(.x$annot.p))*100,2))
purrr::map(tabletemp,~(round(table(.x$annot.p))))

map(tabletemp,~as.data.frame.matrix(table(.x$annot.p,.x$SNP_strand)))%>%
  bind_rows(.id="source")%>% as.data.frame()%>%
  filter(rowSums(across(where(is.numeric)))!=0)%>%
  fwrite(.,file="CLINVAR_annotation_SNV_count.tsv",sep="\t",row.names = TRUE)

map(tabletemp,~as.data.frame.matrix(round(prop.table(table(.x$annot.p,.x$SNP_strand),1)*100,2)))%>%
  bind_rows(.id="source")%>% as.data.frame()%>%
  filter(rowSums(across(where(is.numeric)))!=0)%>%
  fwrite(.,file="CLINVAR_annotation_SNV_prop.tsv",sep="\t",row.names = TRUE)

  #this is again for separate annotation based on gene strand and g4 strand difference
map2(tabletemp,names(tabletemp),~onefunc(as.data.frame(.x),.y))

onefunc<-function(firstlist,DNAorRNA){

  tablelist<-split(firstlist,  firstlist$g4genestrand)

x1<-purrr::map(tablelist,~as.data.frame.matrix(round(prop.table(table(as.character(.x$annot.p),.x$SNP_strand),1)*100,2)))%>% 
  bind_rows(.id="source")%>% as.data.frame()
x2<-purrr::map(tablelist,~as.data.frame.matrix((table(as.character(.x$annot.p),.x$SNP_strand))))%>% 
  bind_rows(.id="source")%>% as.data.frame()

y1<-purrr::map(tablelist,~as.data.frame.matrix(round(prop.table(table(as.character(.x$typemut),.x$which),1)*100,2)))%>% 
  bind_rows(.id="source")%>% as.data.frame()
y2<-purrr::map(tablelist,~as.data.frame.matrix((table(as.character(.x$typemut),.x$which))))%>% 
  bind_rows(.id="source")%>% as.data.frame()

z1<-purrr::map(tablelist,~as.data.frame.matrix(round(prop.table(table(.x$typemut,paste0(.x$V1factor,":",.x$V2factor)),1)*100,2)))%>% 
  bind_rows(.id="source")%>% as.data.frame()
z2<-purrr::map(tablelist,~as.data.frame.matrix(table(.x$typemut,paste0(.x$V1factor,":",.x$V2factor))))%>% 
  bind_rows(.id="source")%>% as.data.frame()

z3<-purrr::map(tablelist,~as.data.frame.matrix(round(prop.table(table(.x$which,paste0(.x$V1factor,":",.x$V2factor)),1)*100,2)))%>% 
  bind_rows(.id="source")%>% as.data.frame()
z4<-purrr::map(tablelist,~as.data.frame.matrix(table(.x$which,paste0(.x$V1factor,":",.x$V2factor))))%>% 
  bind_rows(.id="source")%>% as.data.frame()

fwrite(x2,file=paste0("CLINVAR_annotation_g4genestrand_SNV_count_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)
fwrite(x1,file=paste0("CLINVAR_annotation_g4genestrand_SNV_prop_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)

fwrite(y2,file=paste0("CLINVAR_annotation_typemut_which_count_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)
fwrite(y1,file=paste0("CLINVAR_annotation_typemut_which_prop_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)

fwrite(z1,file=paste0("CLINVAR_typemut3mer_prop_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)
fwrite(z2,file=paste0("CLINVAR_typemut3mer_count_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)
fwrite(z3,file=paste0("CLINVAR_which3mer_prop_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)
fwrite(z4,file=paste0("CLINVAR_which3mer_count_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)

return(list(y1,x1,y2,x2,z1,z2,z3,z4))
}





#======================================================================================================================================================


tiff("CLINVAR_G4_SNV_effect_typmeut_percentage_overall.tiff", units="cm", width=8, height=5, res=300)
cairo_pdf("CLINVAR_G4_SNV_effect_typmeut_percentage_overall.pdf", width=8, height=5, bg = "transparent")

x<-as.data.frame((table(paste0(all_1$typemut),paste0(all_1$V1factor,sprintf('\u2192'),all_1$V2factor))))
  
as.data.frame(prop.table(table(paste0(all_1$typemut),paste0(all_1$V1factor,sprintf('\u2192'),all_1$V2factor)),1)*100) %>% left_join(.,x,by=c("Var1","Var2"))%>%
  filter(Freq.y>0)%>% #   group_by(Var1)%>%
  arrange(desc(Freq.y)) %>%# fwrite("SNV_effect_COSMIC_trinucleotide.context.tsv")
  ggplot(aes(as.factor(Var1), Freq.x,fill=Var2)) +   labs(x="Effect of SNV on 3mer context",y="Frequency",fill="Type of SNV")+
  geom_col( position = "dodge") + 
  geom_text(aes(label=Freq.y,x=as.factor(Var1),y=Freq.x),position = position_dodge(1),vjust=-0.5)+theme_classic()+
  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+theme(legend.text =  element_text(family = "Arial"))

dev.off()

#======================================================================================================================================================



tiff("G4_SNV_effect_which1_percentage.tiff",width= 800, height= 600, units="px", res=300,compression = "lzw")
pdf("G4_SNV_effect_which1_percentage.pdf", width=8, height=6)

x<-as.data.frame((table(paste0(all_1$which),paste0(all_1$V1factor,sprintf('\u2192'),all_1$V2factor))))

as.data.frame(round(prop.table(table(paste0(all_1$which),paste0(all_1$V1factor,sprintf('\u2192'),all_1$V2factor)),1)*100,2))%>% left_join(.,x,by=c("Var1","Var2"))%>%
  filter(Freq.y>0)%>% #   group_by(Var1)%>%
  arrange(desc(Freq.y)) %>%# fwrite("SNV_effect_COSMIC_trinucleotide.context.tsv")
  ggplot(aes(x=as.factor(Var2), y=Freq.x,fill=Var1)) +   labs(x="Effect of SNV on TRINUCLEOTIDE Context",y="Frequency",fill="Type of SNV")+
  geom_col( position = "dodge") + 
  geom_text(aes(label=Freq.y,x=as.factor(Var2),y=Freq.x),position = position_dodge(1),vjust=-0.5)+theme_classic()+
  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+theme(legend.text =  element_text(family = "Arial"))

dev.off()


#--------------------------------------------------------------------------------------------------------------------------------

#======================================================================================================================================================
#======================================================================================================================================================
tabletemp[[1]]

tabletemp<-purrr::map(clinvar_annot.l,~.x%>%
                        dplyr::select(SNP_strand,g4genestrand,annot.p,id,typemut,which,rel_loc,V1factor,V2factor,MFE.x,MFE.y,rel_loc)%>% distinct()%>% 
                        mutate(diff=(MFE.y-MFE.x),g4genestrand=recode(g4genestrand,same="Non-Template", opposite="Template"))%>% as.data.frame())#%>%


library(rstatix)
library(FSA)

library(FSA)
library("dunn.test")
library("DescTools")


quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
stat.test <- map(tabletemp,~.x %>% filter(diff!=0)%>%
  rstatix::group_by(which, g4genestrand) %>% mutate(diff=abs(diff))%>%
  rstatix::wilcox_test(diff ~ annot.p) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()%>% as.data.frame())

stat.test <- map(stat.test,~.x %>% add_xy_position(x = "annot.p", scales = "free"))

stat.test%>% rbindlist(.,idcol="source")%>%as.data.frame()%>% 
  fwrite(.,file="CLINVAR_diffMFE_annotp_strandwhich_wilcox.tsv",sep="\t",quote = F)

plots_MFE_result<-map2(tabletemp,names(tabletemp),~plot_MFE(.x,.y))

#data<-tabletemp[[2]]
#name<-names(tabletemp[[2]])
plot_MFE<-function(data,name){

stat.test <- data %>% filter(diff!=0)%>%
                   rstatix::group_by(which, g4genestrand) %>% mutate(diff=abs(diff))%>%
                   rstatix::wilcox_test(diff ~ annot.p) %>%
                   adjust_pvalue(method = "fdr") %>%
                   add_significance()%>% as.data.frame()

stat.test <- stat.test %>% add_xy_position(x = "annot.p", scales = "free")


p<-ggplot(data = data%>% filter(diff!=0), mapping = aes(reorder(annot.p,-abs(diff)), (diff), fill=annot.p))  +
  guides(fill=F) +
  stat_summary(fun.data = quantiles_95, geom="boxplot")+ facet_grid(which~g4genestrand,drop=TRUE,scales="free")+xlab("Annotation") +
  ylab("ΔMFE ")+
  theme(strip.background = element_rect(fill = "gray90"),
        panel.border = element_rect(colour = "black", fill = NA)) +theme_base( base_family = "Arial")+
  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()  + guides(x = guide_axis(check.overlap = TRUE, angle = 45)) 

#p +
# stat_pvalue_manual(stat.test, hide.ns = TRUE, tip.length = 0.01, inherit.aes = FALSE)+ 
#  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))

p<-p+ stat_compare_means(method = "anova",label = "p.format", method.args = list(test = "Tukey"))
tiff(paste0("CLINVAR_diffMFE_annotp_strandwhich_wilcox",name,".tiff"), units="cm", width=8, height=8, res=300,compression = "lzw")
print(p)
dev.off()

pdf(paste0("CLINVAR_diffMFE_annotp_strandwhich_wilcox",name,".pdf"), width=8, height=8)
print(p)
dev.off()
return(stat.test)
}


#======================================================================================================================================================
#======================================================================================================================================================




base = 6 # set the height of your figure (and font)
expand = 2 # 
j=1
p=list()

p<-purrr::map(tabletemp,~.x %>% 
         dplyr::filter(!is.na(annot.p) )%>% dplyr::select(id,SNP_strand,rel_loc,annot.p,g4genestrand )%>% distinct() %>%
          # filter((SNP_strand %in% c("A->G","G->A","G->T","T->G")))%>% 
         ggplot( aes(x=rel_loc,..count..,fill=factor(g4genestrand))) +
         geom_density( alpha=0.8) +#,fill="#69b3a2", color="#e9ecef",
         scale_x_continuous(breaks = seq(0,1,by=0.3))+theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
         xlab("Relative location of SNV along a G quadruplex")+
         guides(fill=guide_legend("Mutation Type"))+ theme(axis.text.y=element_text(size=rel(1)),axis.text.x=element_text(size=rel(1))) +
         ggtitle("Distribution of SNV across G4")+facet_wrap(~(g4genestrand),scales="free")+
           scale_fill_Publication()+scale_colour_Publication()+theme_Publication())
         
names(p)<-names(tabletemp)
for (i in names(tabletemp)){
  tiff(paste0("CLINVAR_distribution_SNV_g4genestrand",i,".tiff"), units="cm",width=8,height=6,res=300,compression = "lzw")
  print(p[[i]])
  dev.off()
}


#======================================================================================================================================================
#======================================================================================================================================================

swr = function(string, nwrap=20) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

library("ggtext")


map2(tabletemp,names(tabletemp),~annotcount_plot(.x,.y))

#data<-tabletemp[[1]]
#name<-names(tabletemp[[1]])

annotcount_plot<-function(data,name){
tablelist<-split(data,  data$g4genestrand)

myPlots <- list()
for (i in unique(data$SNP_strand[data$SNP_strand %in%  c("A->G","G->A","G->T","T->G","C->G","G->C")])){
  myPlots[[i]] <-  map2(tablelist,names(tablelist),~.x %>% 
                          filter(!is.na(annot.p))%>% filter(SNP_strand==i)%>% group_by(annot.p)%>% filter(n()>5)%>% ungroup()%>%select(id,SNP_strand,rel_loc,annot.p)%>% distinct()%>%# filter(SNP_strand %in% c("A->G","G->A","G->T","T->G") )%>% 
                          ggplot(aes(x=rel_loc,..count..,fill=factor(annot.p))) +
                          geom_density( alpha=0.8) +#,fill="#69b3a2", color="#e9ecef",
                          scale_x_continuous(breaks = seq(0,1,by=0.3))+theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+xlab("Relative location of SNV along a G quadruplex")+
                          guides(fill=guide_legend("Mutation Type"))+  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+
                          ggtitle(swr(paste0(i," ",.y)))+facet_wrap(~annot.p,scales="free_y"))}
#                          theme(plot.title = element_textbox_simple(),axis.text.y=element_text(size=rel(0.2)),axis.text.x=element_text(size=rel(0.2))) }


#cowplot::plot_grid(plotlist =myPlots[[1]],align = "h")

#plot_grid(plotlist=mget(paste0("pl_", 1:10)))

# Create grid
result<-map(1:length(myPlots),~ggpubr::ggarrange(myPlots[[.x]]$`Non-Template`,myPlots[[.x]]$`Template`, # list of plots
                                                 labels = "AUTO", # labels
                                                 common.legend = T, # COMMON LEGEND
                                                 legend = "bottom", # legend position
                                                 align = "h", # Align them both, horizontal and vertical
                                                 ncol = 2))  # number of rows
names(result)<-names(myPlots)

#result1<-result[names(myPlots) %in% c("A->G","G->A","G->T","T->G")]

for (i in 1:length(result)) {
   tiff(paste0("CLINVAR_SNV_g4genestrand_annot_",name,i,".tiff"), units="cm", width=12, height=8, res=300,compression = "lzw")
  #ggsave(paste0("CLINVAR_SNV_g4genestrand_annotation","_",i,".svg"), width=16, height=10)
  print(result[[i]])
  dev.off()}


}





#======================================================================================================================================================
#======================================================================================================================================================

#======================================================================================================================================================
#======================================================================================================================================================
purrr::map2(clinvar_annot.l,names(clinvar_annot.l),~tritemp(.x,.y))

tritemp<-function(data,name){

temp<-data%>% select(id,mref,malt,g4genestrand)%>% distinct()%>%
  mutate(g4genestrand=recode(g4genestrand,same="Non-Template", opposite="Template"))%>%
  as.data.frame()#%>%
object_count<-as.data.frame(table(paste0(temp$mref,":",temp$malt),temp$g4genestrand))%>% 
  filter(Freq>0)%>%    group_by(Var1)%>%
  arrange((Freq))
#fwrite(.,file="tempfile.tsv",sep="\t",row.names = TRUE) 

object_count$Var1<-as.character(object_count$Var1)
object_count$Var2<-as.character(object_count$Var2)

p<-object_count %>% 
  mutate(group1 = if_else(Freq < 20, "Others", Var1))%>%
  group_by(group1,Var2) %>% filter(Freq>0)%>%
  summarize(avg = mean(Freq), count = n()) %>%
  ungroup() %>%
  mutate(group = if_else(group1 == "Others",
                         paste0("Others (n =", count, ")"),
                         group1)) %>%
  mutate(group = forcats::fct_reorder(group, avg)) %>% arrange(desc(count))%>%
  ggplot() + 
  geom_col(aes(group, avg,fill=Var2)) + 
  geom_text(aes(group, avg, label = round(avg, 0)),angle =90) + ylab("count")+ scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+
  facet_wrap(~Var2,scales="free") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))#coord_flip()

tiff(paste0("CLINVAR_tri_g4genestrand_",name,".tiff"), units="cm", width=14, height=6, res=300,compression = "lzw")
print(p)
dev.off()

ggsave(paste0("CLINVAR_tri_g4genestrand_",name,".svg"), width=10, height=6)
print(p)
dev.off()
#object_count$obj_num<- seq(1:nrow(object_count)) 
#object_count %>%
#  mutate(cuml = cumsum(Freq)) %>%
#  ggplot(aes(obj_num)) +
#  geom_tile(aes(y = Freq,# + lag(cuml, default = 0),
 #               height = Freq))


}

#======================================================================================================================================================
plot1<-map(clinvar_annot.l,~.x%>% rowwise%>%
             mutate(g4genestrand=recode(g4genestrand,same="Non-Template", opposite="Template"))%>%
             mutate(gquad.x=as.integer(strsplit(g4_type.x, ':')[[1]][1]),
                    gquad.y=as.integer(strsplit(g4_type.y, ':')[[1]][1])) %>% mutate(across(contains(c("gquad.x","gquad.y")),replace_na,0))%>%
             mutate(gquad=max(gquad.x,gquad.y))%>%as.data.frame%>% filter(!gquad %in% c(1,2,3,11,12,13))) # %>%
  # filter((SNP_strand %in% c('A->G','G->A','G->T','T->G')))%>%


#----------------------------------------------------------------------------------------
# New facet label names for supp variable
#supp.labs <- c("Intergenic", "Non template Strand","Template Strand")
#names(supp.labs) <- c("notapplicable", "opposite","same")
#supp.labs<-swr(supp.labs)

mypalette<-c('#e6194B', '#3cb44b', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#469990', '#9A6324', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#dcbeff', '#fabed4', '#fffac8')
SNVs_c<-unique(plot1[[1]]$SNP_strand)
temp<-map(plot1, ~ .x %>% dplyr::select(id, gquad, SNP_strand, g4genestrand, annot.p)%>%distinct()%>%
            filter(!is.na(g4genestrand) & (annot.p!="CpG Islands"))%>%
            dplyr::select(id,gquad,SNP_strand,g4genestrand)%>% distinct()%>%
            filter(gquad %in% c(4,5,6,7,8)) %>% count(SNP_strand,gquad,g4genestrand) %>%group_by(gquad,g4genestrand) %>%          # now required with changes to dplyr::count()
            mutate(prop = prop.table(n))%>%as.data.frame()%>%
  ggplot( aes(x = gquad, y = prop, color =reorder(SNP_strand,-prop), group = reorder(SNP_strand,-prop))) +
  geom_line() + xlab("Number of Quartets in G quadruplex") + ylab("Proportion of SNV")+ labs(color="SNV in G4")+
  facet_wrap(g4genestrand~.)+
  geom_point()+scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+scale_color_manual(values = mypalette)) 

for (i in 1:length(temp)){
tiff(paste0("CLINVAR_G4_percentage_quartet_linegraph_g4genestrand_",names(temp)[[i]],".tiff"), units="cm", width=8, height=5, 
     res=300,compression = "lzw")
print(temp[[i]])
dev.off()
pdf(paste0("CLINVAR_G4_percentage_quartet_linegraph_g4genestrand_",names(temp)[[i]],".pdf"), width=8, height=5)
print(temp[[i]])
dev.off()
  
}





#G4_expr$gquad<-as.factor((as.integer(G4_expr$gquad)))
library(forcats)

temp1<-map(plot1,~.x%>%dplyr::select(id,gquad,SNP_strand,g4genestrand,annot.p)%>%distinct()%>%filter(!is.na(g4genestrand) & (annot.p!="CpG Islands"))%>%
  dplyr::select(id,gquad,SNP_strand,g4genestrand)%>% distinct()%>%
  filter(gquad %in% c(4,5,6,7,8)) %>% count(SNP_strand,gquad,g4genestrand) %>%group_by(gquad,g4genestrand) %>%          # now required with changes to dplyr::count()
  mutate(prop = prop.table(n))%>%as.data.frame())

plot2<-temp1[[1]]%>%#arrange(desc(n))%>%
  ggplot(., aes( y=n, x=gquad,fill=reorder(SNP_strand,-prop),  group = reorder(SNP_strand,-prop))) + 
  geom_bar(position="stack", stat="identity")+ xlab("Number of Quartets in G quadruplex") + ylab("Counts of SNV")+ labs(color="SNV in G4")+
  facet_wrap(~g4genestrand)+
  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+
  guides(fill = guide_legend(title = "SNV in G quadruplex"))



plot2


library(ggpubr)
tiff("Clinvar_G4_percentage_quartet_linegraph_barplot_g4genestrand.tiff", units="cm", width=15, height=15, res=300)
ggsave("Clinvar_G4_percentage_quartet_linegraph_barplot_g4genestrand.pdf", width=8, height=10)

cowplot::plot_grid(temp[[1]],plot2,ncol=1,align="v")

#ggarrange(plot1, plot2,
#          align='v', ncol = 2,labels=c('A', 'B'),
#          common.legend = T)
dev.off()


#======================================================================================================================================================
map2(tabletemp,names(tabletemp),~plot2(.x,.y))

plot2<-function(data,name){
  tablelist<-split(data,  data$g4genestrand)
  
  
p<-map(tablelist,~as.data.frame(round(prop.table(table(.x$annot.p,.x$which),1)*100,2)))%>% 
  bind_rows(.id="source")%>% as.data.frame()%>% 
  filter(rowSums(across(where(is.numeric)))!=0)%>% mutate(Effect=Var2)%>%mutate(Var=NULL)%>%  filter(rowSums(across(where(is.numeric)))!=0) %>% drop_na(Var1) %>%
  ggplot( aes(x = Var1, y = Freq, fill = Effect)) + facet_grid(~source, drop=TRUE, scales="free"  )+
  geom_bar(position="stack", stat="identity")+  xlab("SNV in G quadruplex region across CLINVAR database") + ylab("Percentage of SNV")+
  labs(color="SNV in G4")+  theme_classic( base_family = "Arial") +
  scale_fill_Publication()+scale_colour_Publication()+theme_Publication() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

tiff(paste0("Clinvar_G4_prop_barplot_effect",name,".tiff"), units="cm", width=8, height=5, res=300)
print(p)
dev.off()
pdf(paste0("Clinvar_G4_prop_barplot_effect",name,".pdf"), width=8, height=5)
print(p)
dev.off()
return(p)
}


map2(tabletemp,names(tabletemp),~plot3(.x,.y))


plot3<-function(data,name){
  tablelist<-split(data,  data$g4genestrand)
  
  
  p<-map(tablelist,~as.data.frame(table(.x$annot.p,.x$which)))%>%
    bind_rows(.id="source")%>% as.data.frame()%>% 
    dplyr::filter(rowSums(across(where(is.numeric)))!=0)%>% mutate(Effect=Var2)%>%mutate(Var=NULL)%>%  filter(rowSums(across(where(is.numeric)))!=0) %>% drop_na(Var1) %>%
    ggplot( aes(x = Var1, y = Freq, fill = Effect)) + facet_grid(~source, drop=TRUE, scales="free"  )+
    geom_bar(position="stack", stat="identity")+  xlab("Annotation") + ylab("Count of SNV")+
    labs(color="SNV in G4")+theme_base(base_family = "Arial")+
    scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+coord_flip()
  tiff(paste0("Clinvar_G4_count_barplot_effect",name,".tiff"), units="cm", width=8, height=5, res=300)
  print(p)
  dev.off()
  pdf(paste0("Clinvar_G4_count_barplot_effect",name,".pdf"), width=8, height=5)
  print(p)
  dev.off()
  return(p)
}




#======================================================================================================================================================



library("MKinfer")
combn(levels(tabletemp[[1]]$annot.p), 2, function(x) {
  perm.t.test(diff~annot.p,data = subset(tabletemp[[1]], annot.p %in% x), nperm=100, paired = F)
}, simplify = FALSE) -> result

result

#======================================================================================================================================================
#======================================================================================================================================================











#======================================================================================================================================================
#======================================================================================================================================================

#rm(G4_var_GSM,bedA,GSM_hg38,result,GSM,chainObject,x,x1)
#======================================================================================================================================================
#check if any variants is repeated



x2%>% group_by(name)%>%
  summarise(n=n())%>% arrange(desc(n))%>%head(7)
#if not the below lines not required
x<-x2%>% group_by(name)%>%
  summarise(n=n())%>% arrange(desc(n))%>%head(7)%>% select(name)%>% unlist(.,use.names = F)#as.data.frame()
x2[x2$name %in% x,]


#======================================================================================================================================================
#======================================================================================================================================================
library(annotatr)
annots <- builtin_annotations()[grep("hg38",builtin_annotations())]

annots=c("hg38_genes_1to5kb"  , "hg38_genes_promoters"  ,               
         "hg38_genes_5UTRs",  "hg38_genes_exons","hg38_genes_firstexons" ,         
         "hg38_genes_introns"   , "hg38_genes_cds" ,"hg38_genes_intronexonboundaries" ,"hg38_genes_exonintronboundaries",
         "hg38_genes_3UTRs","hg38_genes_intergenic","hg38_cpg_islands"   ,   
         "hg38_enhancers_fantom", "hg38_lncrna_gencode", "hg38_basicgenes" )

# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg38',annots)


seqlevels(annotations, pruning.mode="coarse") <- paste0('chr', seq(1,22))

G4_gm<-makeGRangesFromDataFrame(x2,seqnames.field = "left",start.field = "right",
                                end.field = "right",strand.field = "G_quadstrand",ignore.strand = F,keep.extra.columns = TRUE)

# Randomize the input regions
genome<-'hg38'

genome(G4_gm)<-'hg38'
#Intersect the regions we read in with the annotations
dm_annotated = annotate_regions(
  regions = G4_gm,
  annotations = annotations,
  ignore.strand = T,
  quiet = FALSE)
# A GRanges object is returned
print(dm_annotated)
genome(dm_annotated)<-'hg38'

# Coerce to a data.frame
df_dm_annotated = data.frame(dm_annotated)

# See the GRanges column of dm_annotaed expanded
print(head(df_dm_annotated))

genome(G4_gm) <- "hg38"

# Randomize the input regions
dm_random_regions = randomize_regions(
  regions = G4_gm,
  allow.overlaps = T,
  per.chromosome = T)

# Annotate the random regions using the same annotations as above
# These will be used in later functions
dm_random_annotated = annotate_regions(
  regions = dm_random_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = TRUE)


#======================================================================================================================================================
#======================================================================================================================================================
#----------------------------------------------------------------------------------------------------------------------------------

clinvar_annotated<-load("C:/Users/a0neup03/OneDrive - University of Louisville/Documents/SNV_paper/g4snv/clinvar_g4_annotated_usethis.RData")

head(df)
df$g4genestrand<-ifelse (df$annot.strand==df$strand,"same","opposite")
df$g4genestrand<-ifelse (df$annot.strand=="*","notapplicable",df$g4genestrand)
#clinvar_annotated[as.numeric(G4_expr$right) %in% as.numeric(clinvar_annotated$start)  ,]


nrow(df)

#df[G4_expr$id %in% df$name,]%>% select(name)%>% distinct()%>% nrow()

G4_expr<-left_join(G4_expr,all_1[,c("id","V1factor","V2factor","mref","malt")],by="id")
clinvar_annot<-left_join(G4_expr,df[,c("name","annot.gene_id", "annot.symbol","annot.type","annot.p","g4genestrand")],by=c("id"="name"),multiple="all")%>%  as.data.frame()


clinvar_annot$which.ED<-ifelse(clinvar_annot$ED.x > clinvar_annot$ED.y,"ref_more","alt_more")
clinvar_annot$which.ED<-ifelse(clinvar_annot$ED.x==clinvar_annot$ED.y,"no change",clinvar_annot$which.ED)
transitions<-c("A->G","G->A","C->T","T->C")
transversions<-c("A->T","A->C","C->A","C->G","T->A","T->G","G->C","G->T")

clinvar_annot$typemut<-ifelse(clinvar_annot$SNP %in% transitions,"Transition","Transversion")



#now we separate RNA and DNA annotations
#clinvar_annot<-


  clinvar_annot_a<-clinvar_annot %>% filter(annot.type %in% c("hg38_genes_promoters", "hg38_genes_intergenic")) %>% 
    # Keep only the rows with specific value
    dplyr::slice(rep(1:n(), each = 2)) %>% # Repeat each row twice
    mutate(annot.group = ifelse(annot.type == "hg38_genes_promoters", c("DNA","RNA"), c("DNA","RNA")))%>%as.data.frame()
  
  
  clinvar_annot_b<- clinvar_annot%>% filter(!annot.type %in% c("hg38_genes_promoters", "hg38_genes_intergenic")) %>% 
    mutate(annot.group=case_when((annot.type=="hg38_lncrna_gencode")~"DNA",
                         (annot.type=="hg38_enhancers_fantom")~"DNA",
                         (annot.type=="hg38_genes_1to5kb")~"DNA",
                         (annot.type=="hg38_cpg_islands")~"DNA",
                         (annot.type=="hg38_genes_introns")~"DNA",
                         (annot.type=="hg38_genes_cds")~"RNA",
                         (annot.type=="hg38_genes_exons")~"DNA",
                         (annot.type=="hg38_genes_5UTRs")~"RNA",
                         (annot.type=="hg38_genes_3UTRs")~"RNA",
                         (annot.type=="hg38_genes_exonintronboundaries")~"DNA",
                         (annot.type=="hg38_genes_intronexonboundaries")~"DNA"))

clinvar_annot<-  rbind(clinvar_annot_b,clinvar_annot_a)
clinvar_annot<-clinvar_annot %>% arrange(desc(seqnames),desc(start.x),annot.gene_id) %>% group_by_at(vars(-annot.gene_id,-annot.symbol,-annot.type,-annot.p))%>% tidyr::fill(annot.gene_id,annot.symbol)%>% ungroup()

clinvar_annot<-setDT(clinvar_annot)

clinvar_annot.l<-split(clinvar_annot,clinvar_annot$annot.group)


clinvar_annot.l<-purrr::map2(clinvar_annot.l,names(clinvar_annot.l),~get_1annot(.x,.y))
(table(clinvar_annot.l[[2]]$annot.p1))
(table(clinvar_annot.l[[1]]$annot.p1))

table(clinvar_annot.l[[2]]$n)

clinvar_annot.l[[1]]%>% arrange(desc(n))

get_1annot<-function(dat1,name){

  if (name=="DNA"){
    dat1$annot.p1<- factor(dat1$annot.p,  levels = c("EXON", "INTRON", "PROMOTERS",
                                                              "Enhancers", "CpG Islands",
                                            "lncRNA GENCODE","Nearby1to5kb","Intergenic","EXON/INTRON Boundaries","INTRON/EXON Boundaries", "hg38_genes_firstexons" ))
  }  else if (name=="RNA"){
    dat1$annot.p1<- factor(dat1$annot.p,  levels = c("CDS", "5' UTR", "3' UTR","PROMOTERS","Intergenic"))
  }
  
dat1<-dat1%>% arrange(annot.p1)%>% as.data.frame()%>%group_by_at(vars(-annot.gene_id,-annot.symbol))%>%
  summarise(genes.id=paste0(unique(annot.gene_id[!is.na(annot.gene_id)]),collapse=", "),
         genes.symbol=paste0(unique(annot.symbol[!is.na(annot.symbol)]),collapse=", "),
         n=n_distinct(annot.gene_id,na.rm=TRUE))%>%
  arrange(desc(n))%>% as.data.frame()
#Should retain the first unique row for each id
dat1<-dat1[order(dat1$id, dat1$annot.p1,-dat1$n),]

dat1<-dat1[!duplicated(dat1$id),]

return(dat1)
}




dat1<-clinvar_annot.l[[1]]


library(rstatix)
library(FSA)
dat1$MFediff<-(dat1$MFE.y-dat1$MFE.x)

res.aov<-dat1       %>% filter(annot.p!="Nearby1to5kb")%>% 
  filter(MFediff!=0)%>%select(id,annot.p,MFediff,which,g4genestrand,typemut,strand)%>%distinct()%>% group_by(which,g4genestrand)%>%
  filter(!is.na(annot.p)) %>% anova_test(  MFediff ~ annot.p,wid=id)
res.aov
res.aov<-dat1       %>% filter(annot.p!="Nearby1to5kb")%>% 
  filter(MFediff!=0)%>%select(id,annot.p,MFediff,which,g4genestrand,typemut,strand)%>%distinct()%>% group_by(typemut,g4genestrand)%>%
  filter(!is.na(annot.p)) %>% anova_test(  MFediff ~ annot.p,wid=id)
res.aov


dat2<-
dat1 %>% filter(annot.p1!="Nearby1to5kb")%>% 
  select(id,annot.p,MFediff,which,MFE.x,MFE.y,g4genestrand,strand,typemut) %>% distinct() %>% filter(!is.na(annot.p))

wilcox.test(dat2$MFE.x[(dat2$which=="Further stabilized")],dat2$MFE.y[(dat2$which=="Further stabilized")])
t.test(dat2$MFE.x[(dat2$which=="Further stabilized")],dat2$MFE.y[(dat2$which=="Further stabilized")])
wilcox.test(dat2$MFE.x[(dat2$which=="Destabilized")],dat2$MFE.y[(dat2$which=="Destabilized")])
t.test(dat2$MFE.x[(dat2$which=="Destabilized")],dat2$MFE.y[(dat2$which=="Destabilized")])


res.posthoc <- dat2 %>% filter(annot.p!="Nearby1to5kb")%>% filter(which!="no change")%>%
   select(id,annot.p,MFediff,which,g4genestrand) %>% distinct() %>% filter(!is.na(annot.p)) %>% group_by(which)%>%
  pairwise_wilcox_test(MFediff ~ annot.p, p.adjust.method = "bonferroni", paired = F, ref.group = "PROMOTERS")

res.posthoc 

dat1 %>% filter(annot.p1!="Nearby1to5kb")%>% 
  select(id,annot.p,MFediff,which,g4genestrand) %>% distinct() %>% group_by(g4genestrand,which)%>%filter(!is.na(annot.p)) %>% kruskal.test(MFediff ~ annot.p)

#dat2<-

  dat2 %>% filter(annot.p!="Nearby1to5kb")%>% #filter(which=="Further stabilized")%>%
  select(id,annot.p,MFediff,which,g4genestrand) %>% distinct() %>% filter(!is.na(annot.p))%>% group_by(g4genestrand)%>%
    dunn_test(MFediff~ annot.p,p.adjust.method = "bonferroni")%>% as.data.frame()#%>% fwrite(.,file="CLINVAR_dunntest_genestrand_annot_overall.tsv",sep='\t',col.names=T,quote=F)
  
  
  dat2 %>% filter(annot.p!="Nearby1to5kb")%>% filter(which=="Destabilized")%>%
    select(id,annot.p,MFediff,which,strand,typemut) %>% distinct() %>% filter(!is.na(annot.p))%>% group_by(typemut,strand)%>%
    dunn_test(MFediff~ annot.p,p.adjust.method = "bonferroni")%>% as.data.frame()#%>% fwrite(.,file="CLINVAR_dunntest_genestrand_annot_overall.tsv",sep='\t',col.names=T,quote=F)
  
  
#dunn test
dunn.test(x=abs(dat2$MFediff[(dat2$which=="Further stabilized")]),dat2$annot.p[dat1$which=="Further stabilized"],method = "bonferroni")
dunn.test(dat2$MFediff[(dat2$which=="Destabilized")],dat2$annot.p[dat2$which=="Destabilized"],method = "bonferroni")

dat2 %>% filter(annot.p!="Nearby1to5kb")%>% filter(which=="Destabilized")%>% 
  select(id,annot.p,MFediff,which,g4genestrand) %>% distinct() %>% filter(!is.na(annot.p))%>%group_by(g4genestrand)%>% 
  #kruskal_test(MFediff~ annot.p)
  dunn_test(MFediff~ annot.p,p.adjust.method = "bonferroni")

dat2 %>% filter(annot.p!="Nearby1to5kb")%>%# filter(which=="Further stabilized")%>%
  group_by(g4genestrand)%>% 
  select(id,annot.p,MFediff,which,g4genestrand) %>% distinct() %>% filter(!is.na(annot.p))%>%
  dunn_test(MFediff~ annot.p,p.adjust.method = "BY")

#dat2<-
#  dat1 %>% filter(annot.p1!="Nearby1to5kb")%>% 
# select(id,annot.p,MFediff,which,MFE.x,MFE.y,G_quadstrand) %>% distinct() %>% filter(!is.na(annot.p))

dat2%>% 
   distinct() %>% filter(!is.na(annot.p))%>% filter(which=="Destabilized")%>%
  select(id,annot.p,MFediff,which,g4genestrand) %>% distinct() %>% filter(!is.na(annot.p))%>% group_by(g4genestrand)%>%
  dunn_test(MFediff~ annot.p,p.adjust.method = "bonferroni")#%>% as.data.frame()%>% fwrite(.,file="CLINVAR_dunntest_genestrand_annot_destabilizd.tsv",sep='\t',col.names=T,quote=F)

dat1%>%  filter(annot.p!="Nearby1to5kb")%>% filter(which=="Destabilized")%>%
  select(id,annot.p,MFediff,which,strand) %>% distinct() %>%   filter(!is.na(annot.p))%>% group_by(strand)%>%
  dunn_test(MFediff~ annot.p,p.adjust.method = "bonferroni")#%>% as.data.frame()%>% fwrite(.,file="CLINVAR_dunntest_G4strand_annot_destabilizd.tsv",sep='\t',col.names=T,quote=F)


ggboxplot(dat1, x = "annot.p1", y = "MFE.x", color = "annot.p1")
library(viridis)
# Represent it
p <- dat1 %>% filter(which!="no change")%>%
  ggplot( aes(x=log(abs(MFediff)), fill=which)) +
  geom_histogram( color="#e9ecef", alpha=1, position = 'dodge') + scale_fill_viridis(discrete=TRUE) +
  labs(fill="")
p

table(dat1$annot.p1,dat1$genes.id)
dat1%>% mutate(MFEabs=log(abs(MFediff)))%>%
  pairwise_t_test(MFEabs ~ annot.type, p.adjust.method = "bonferroni")%>% filter(p.adj<0.05)







library(Qtlizer)

qt<-G4_expr%>% dplyr::select(seqnames,start.x,end.x,width,strand,id)%>%distinct()%>% 
  reframe(seqnames=paste0("hg38:",gsub("chr","",seqnames),":",start.x))#%>%fwrite(.,"cosmic_qtlizer.input.txt",col.names = F,row.names = F,quote = F,sep="\t")
qt$seqnames
te_qt<-split(qt$seqnames,             # Applying split() function
             ceiling(seq_along(qt$seqnames) / 48))

qtlizer_results<-purrr::map(te_qt,~testit(.x,times=1))

testit <- function(x,times){
  Sys.sleep(times)
  granges = get_qtls(x, return_obj = "granges", ref_version = "hg38",max_terms = 1)
  Sys.sleep(times)
  times<-ifelse(times>50,10,times)
  times=times+10
  return(as.data.frame(granges))
}


df = get_qtls(qt$seqnames, corr = 0.8, ld_method = "r2")
granges


clinvar_annot.l





library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)     # hg38 genome
library(BSgenome)


available.SNPs()
library(bedr)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)


snps.bed.file<-G4_expr%>% dplyr::select(seqnames,right,right,width,strand,id)%>%distinct()%>% 
  reframe(seqnames=as.character(gsub("chr","",seqnames)),
          start=as.numeric(right),
          end=as.numeric(right))%>% reframe(name=paste0(seqnames,":",start,"-",end))%>% fwrite(.,file="CLINVAR_kaviar.tsv",sep="\n",row.names = F,col.names = F,quote = F)

clinVar.vcf.example<-read.vcf("clinvar_kaviar_annotated.txt")
clinvar_rsid<-unique(clinVar.vcf.example$vcf$ID)

system.time(motifR_result<-lapply(clinvar_rsid,function(x)motifRtry(x)))

motifRtry<-function(clinvarrsid){
snps.mb <- snps.from.rsid(rsid = clinvar_rsid[!is.na(clinvar_rsid)],
                          dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38)


results1 <- motifbreakR(snpList = snps.mb, filterp = T,
                                    pwmList = motifbreakR_motif,
                                    threshold = 0.01,
                                    method = "default",
                                    bkg = c(A=0.15, C=0.2, G=0.5, T=0.15),
                                    BPPARAM =BiocParallel::SerialParam())

#rs1006140 <- results[names(results) %in% "rs1006140"]
result<-as.data.frame(results1@elementMetadata)
return(result)
}

motifbreakr.results <- motifbreakR(snpList = variants, pwmList = MotifDb, threshold = 0.9)
plotMB(results = motifbreakr.results, rsid = "rs7837328", effect = "strong")


library(MotifDb)
table(mcols(MotifDb)$organism, mcols(MotifDb)$dataSource)

data(motifbreakR_motif)

library(BiocParallel)
multicoreParam <- SnowParam(workers = 8)
# process in parallel

system.time(results1 <- motifbreakR(snpList = snps.mb, filterp = T,
                                    pwmList = motifbreakR_motif,
                                    threshold = 0.01,
                                    method = "default",
                                    bkg = c(A=0.15, C=0.2, G=0.5, T=0.15),
                                    BPPARAM =BiocParallel::SerialParam()))

#rs1006140 <- results[names(results) %in% "rs1006140"]
result<-as.data.frame(results1@elementMetadata)
result%>% group_by(SNP_id,effect)%>% summarise(n=n(),n_distinct(geneSymbol),n_distinct(providerId),adiff=max(abs(alleleEffectSize)))%>% arrange(desc(adiff))%>% head()

exm_mut<-"chr10:122143482:G:A"
exm_mut <- results1[names(results1) %in% exm_mut]

exm_mut<-results1[(results1$alleleEffectSize) > 0.2]
exm_mut_p <- calculatePvalue(exm_mut)
tiff("G4_mut_chr10_122143482_G_A.tiff", units="in", width=8, height=5, res=300)
a<-plotMB(exm_mut, rsid = "chr10:122143482:G:A", effect = "strong", reverseMotif = TRUE)
dev.off()


exm_mut<-"chr22:21851158:G:A"
exm_mut <- results1[names(results1) %in% exm_mut]

exm_mut_p <- calculatePvalue(exm_mut)
exm_mut<-results1[(results1$alleleEffectSize) > 0.25]

tiff("G4_mut_chr22_21851158_G_A.tiff", units="in", width=8, height=5, res=300)
plotMB(exm_mut, rsid = "chr22:21851158:G:A", effect = "strong", reverseMotif = TRUE)
dev.off()

chr19:5587257:T:G
chr18:75285891:A:C
chr11:191891:T:A
chr10:122143482:G:A
chr1:201665601:T:C
library(rtracklayer)
snps <- rtracklayer::import(file=snps.bed.file, format = "bed")
nrow(as.data.frame(results1))
#http://gong_lab.hzau.edu.cn/PancanQTL/download
G4_expr%>% filter(start.x==1295113)

