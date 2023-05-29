
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


cb <- function(df, sep="\t", dec=",", max.size=(200*1000)){
    # Copy a data.frame to clipboard
    fwrite(df, paste0("clipboard-", formatC(max.size, format="f", digits=0)), sep=sep, row.names=FALSE, dec=dec)
  }




clinvarseq<-readRDS("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/clinvar.rds")
clinvar<-readRDS("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/clinvar_1.RDS")

clinvar1<-left_join(clinvar,clinvarseq[,c("G4_id","REF","ALT","loc_start")],by="G4_id")


clinvar_stats<-clinvar[, .(G4hunter_ref = mean(V9), G4hunter_alt = mean(V13),count=.N,sd=sd(V9),sd2=sd(V13)), by = .(V12,V16)]%>% fwrite(file="CLINVAR_g4hunter.tsv")


clinvar1<-clinvar1%>% filter(V12==1 | V12==-1 | V16==1 | V16== -1)
clinvar1<-setDT(clinvar1)
clinvar1[, G4_refrc := ifelse((V12 == -1 | V16 == -1) , sapply( G4_ref, function(x) as.character(
reverseComplement(DNAString(x)))), sapply(G4_ref,function(x) as.character(x)))]


clinvar1[, G4_altrc:= ifelse(V12 == -1 | V16 == -1, sapply( G4_alt, function(x) as.character(
reverseComplement(DNAString(x)))), sapply(G4_alt,function(x) as.character(x) ))]
clinvar1$nameimp<-clinvar1$name

clinvar1$name<-gsub("_G4_id","",clinvar1$name)
clinvar1<-separate(data = clinvar1, col = name, into = c("left", "right"), sep = "_")

 clinvar1$start<-as.numeric(clinvar1$right)-30
clinvar1$end<-as.numeric(clinvar1$right)+30

clinvar1$name<-paste0(clinvar1$left,":",clinvar1$start,"-",clinvar1$end,"||","loc:",clinvar1$loc_start,"||",clinvar1$REF,"&",clinvar1$ALT)
table(clinvar1$V12,clinvar1$V16)


#writeFastaref(clinvar1,"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/clinvar_1_ref_g4.fasta")
#writeFastaalt(clinvar1,"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/clinvar_1_alt_g4.fasta")



ref_RNAFOLD_clinvar<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/RNAfold/clinvar_1_ref_g4/clinvar_1_ref_g4RNAFOLD_allG4COSMICfasta_var_features.txt")
alt_RNAFOLD_clinvar<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/RNAfold/clinvar_1_alt_g4/clinvar_1_alt_g4RNAFOLD_allG4COSMICfasta_var_features.txt")
RNAFOLD_clinvar<-left_join(ref_RNAFOLD_clinvar[,c("MFE","MFEadj","EFE","CDE","MEAFE","MEA","ED","MFEadj.GC","MFEadj.BP","MEAFEadj","EDadj","CFE","EFEadj","name")],alt_RNAFOLD_clinvar[,c("MFE","MFEadj","EFE","CDE","MEAFE","MEA","ED","MFEadj.GC","MFEadj.BP","MEAFEadj","EDadj","CFE","EFEadj","name")],by="name")






ref_RNAFOLD<-setDT(ref_RNAFOLD_clinvar)
alt_RNAFOLD<-setDT(alt_RNAFOLD_clinvar)
ref_RNAFOLD<- unique(ref_RNAFOLD, by = c('name'))

alt_RNAFOLD<- unique(alt_RNAFOLD, by = c('name'))
RNAFOLD<-merge(ref_RNAFOLD[,c("ID","MFE","MFEadj","EFE","CDE","MEAFE","MEA","ED","MFEadj.GC","MFEadj.BP","MEAFEadj","EDadj","CFE","EFEadj","name")],alt_RNAFOLD[,c("MFE","MFEadj","EFE","CDE","MEAFE","MEA","ED","MFEadj.GC","MFEadj.BP","MEAFEadj","EDadj","CFE","EFEadj","name")],by="name")
head(RNAFOLD)



#
whichdiff<-function(df,col){
colref<-paste0(col,".x")
colalt<-paste0(col,".y")
colnamedf<<-c(colref,colalt)
print(colnamedf)
df<-as.data.frame(df[,colnamedf])
df[,"which"]<-ifelse(df[,colref] > df[,colalt],"alt_sml","ref_sml")
df[,"which"]<-ifelse(df[,colref] ==df[,colalt],"no change",df[,"which"])
print(table(df$which))
return(cbind(RNAFOLD,df[!names(df) %in% colnamedf]))
}

df1<-whichdiff(as.data.frame(RNAFOLD),"MFE")

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

levels(RNAFOLD$which)<-c("Further stabilized","no change","Destabilized")
#for ED
#levels(RNAFOLD$which)<-c("less diversity","no change","more diversity")

#RNAFOLD$ID<-gsub('(?<=[||]|[||])\\w+', '', RNAFOLD$ID, perl = TRUE)
RNAFOLD$ID<-gsub("[||].*[||]", "||", RNAFOLD$ID)


RNAFOLD<-separate(
  RNAFOLD,
  col=name,
  into=c("region","loc","SNP"),
  sep = "\\|\\|",
  remove = FALSE)
head(RNAFOLD)

RNAFOLD$SNP<-gsub("&","->",RNAFOLD$SNP)
clinvar1$ID<-gsub("[||].*[||]", "||", clinvar1$name, perl = TRUE)

merged<-merge(clinvar1,RNAFOLD,by=c("ID","name"))
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


data_df<-test# %>% filter(nchar(SNP)==4)



writeFastaref(test,"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/clinvar_1_ref_g4.fasta")
writeFastaalt(test,"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/clinvar_1_alt_g4.fasta")



temp_clinvar_r<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/quadparser_clinvar_ref_g4.tsv",sep="\t")
temp_clinvar_a<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/quadparser_clinvar_alt_g4.tsv",sep="\t")

colnames(temp_clinvar_r)[c(2,3,7,5)]<-c("g4_start","g4_stop","g4seq","g4_type")
colnames(temp_clinvar_a)[c(2,3,7,5)]<-c("g4_start","g4_stop","g4seq","g4_type")


temp_clinvar_a%>% group_by(V1)%>% mutate(n=n())%>% arrange(desc(n))%>% head()
temp_clinvar_a<-temp_clinvar_a[!duplicated(temp_clinvar_a), ]
temp_clinvar_r<-temp_clinvar_r[!duplicated(temp_clinvar_r),]
rm(temp_coding_r,temp_coding_a,temp_NCV_a,temp_NCV_r)




data1<-left_join(data_df,temp_clinvar_r[,c("V1","g4_start","g4_stop","g4seq","g4_type")],by=c("name"="V1"))
data1na<-data1%>% filter(is.na(g4_start))%>% select(-c(g4_start,g4_stop,g4seq,g4_type))
data1<-data1%>% filter(!is.na(g4_start))

data1na<-left_join(data1na,temp_clinvar_a[,c("V1","g4_start","g4_stop","g4seq","g4_type")],by=c("name"="V1"))

data_total<-rbind(data1,data1na)
rm(data1  ,data1na)

#data_total$affectg4<-ifelse(data_total$g4_start<31 &data_total$g4_stop>29,"affectsg4","doesnot")

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


#save(data_total_1,file="/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/full_clinvar_G4_duplicatesremoved.Rdata")
#load("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/full_clinvar_G4_duplicatesremoved.Rdata")
				(affectg4=="affect_g4")~ (30-g4_start)/(g4_stop-g4_start)))
#save(data_total_1,file="/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/full_codingNCV_G4_duplicatesremoved.Rdata")
#load("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/full_codingNCV_G4_duplicatesremoved.Rdata")

#data_total_1%>% group_by(name,G4_refrc,G4_altrc,start,g4seq,rel_loc)%>% mutate(numberofg4=n(),affectg41=paste0(affectg4,collapse="|"),g4seq11=paste0(g4seq,collapse="|"))%>% arrange(desc(numberofg4))%>%head()%>% as.data.frame() 
#data_total_1<-data_total_1[!duplicated(data_total_1$name),]

#data_total_1%>% group_by(name,G4_refrc,G4_altrc,start,end)%>% mutate(n=n())%>% arrange(desc(n))%>%head()
data_total_2<-data_total_1%>% filter(affectg4=="affect_g4")
#temp<-split(data_total_1,data_total_1$affectg4)
temp<-split(data_total_2,data_total_2$G_quadstrand)
temp<-split(data_total_2,data_total_1$affectg4)



#--------------------------------------------------
data_total_2<-setDT(data_total_2)
temp<-split(data_total_2,data_total_2$G_quadstrand)

test_stats<-data_total_2[, .(G4hunter_ref = mean(V9), G4hunter_alt = mean(V13),count=.N,sd=sd(V9),sd2=sd(V13)), by = .(V12,V16)]%>% fwrite(file="CLINVAR_g4hunter.tsv")
data_total_2%>% mutate(chr=gsub(":.*$", "", name),start=start+g4_start,end=start+g4_start+g4_stop,score=".")%>%
select(chr,start,end,name,score,G_quadstrand)%>%
fwrite(.,"G4_aroundSNVCLINVAR_all.bed",sep="\t",col.names=FALSE,quote=FALSE)



data_dfannot<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/G4_aroundSNVCLINVAR_all_annotated.bed",skip=1)



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


temp<-split(data_dfannot1,data_dfannot1$G_quadstrand)
#temp<-split(data_dfannot1,data_dfannot1$V8)

mat<-map(temp,~xtabs(~V8+paste0(real_mut,"::"),.x))

mat<-map(temp,~xtabs(~which1+paste0(real_mut),as.matrix(.x)))

temp<-map(temp,~setDT(.x))
map(temp,~.x[ , .N, by = .(which1, real_mut) ]%>% 
dcast(which1~real_mut,value.var="N"))%>%rbindlist(id="source")%>% cb(.)


map(temp,~.x[ , .N, by = .(V8,which1) ]%>% 
dcast(V8~which1,value.var="N"))%>%rbindlist(id="source")%>%cb(.)
map(temp,~.x[ , .N, by = .(left, real_mut) ]%>% 
dcast(paste0(real_mut,"::")~left,value.var="N"))%>%rbindlist(id="source")%>% fwrite(.,file="clinvar_bychrreal_mut.tsv")

temp_prom<- map(temp,~.x%>% filter(V8=="promoter-TSS"))


map(temp,~.x[ , .N, by = .(V8, real_mut,which1) ]%>% dcast(paste0(real_mut,"::",V8)~which1,value.var="N")#%>% cb(.)
map(temp,~.x[ , .N, by = .(real_mut) ]%>% dcast(paste0(real_mut,"::")~which1,value.var="N")%>% arrange(desc(Destabilized))) #cb(.)
map(temp_prom,~.x[ , .N, by = .(left, real_mut) ]%>% dcast(left~paste0(real_mut,"::"),value.var="N"))#%>% cb(.)

map(temp,~ggplot(.x, aes(x = factor(real_mut), y = rel_loc)) +
  geom_violin() +  geom_jitter(height = 0, width = 0.08, aes(colour = factor(real_mut))) )
  #scale_colour_manual(name="colour", values=c("grey", "purple", "pink")))
 # geom_boxplot(width=0.3))
dev.off()
#(real_mut %in% c("C->A","C->T","A->C") & (G_quadstrand=="-")
 (G_quadstrand=="+")) |

map(temp,~.x %>% # filter((SNP_strand %in% c("G->A","G->T","T->G")))%>%
ggplot(.) +
geom_density(aes(x=rel_loc,y=..count..,fill=SNP_strand), alpha=0.7) + 
  facet_wrap( SNP_strand~.,scales="free")+ theme_bw() +
scale_colour_Publication()+
 labs(x="Relative Location of SNV along a G quadruplex")+
#theme_Publication()+
theme(axis.text=element_text(size=rel(0.7)),axis.title=element_text(size=10),
  panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
       panel.background = element_blank()))

dev.off()



# Make the histogram
map(temp,~.x %>% filter(!is.na(V8))%>%
  ggplot( aes(x=rel_loc,fill=which1)) +
    geom_density(aes(y=..count..), alpha=0.8) +#,fill="#69b3a2", color="#e9ecef",
scale_x_continuous(breaks = seq(0,1,by=0.1))+
    ggtitle("Distribution of SNV across G4")) #+




map(temp,~.x %>% filter(!is.na(V8))%>% filter(V8 %in% c("3' UTR","TTS","exon","promoter-TSS") &
real_mut %in% c("A->G","G->A","C->T","T->C"))%>% 
  ggplot( aes(x=rel_loc,..count..,fill=real_mut)) +
    geom_density( alpha=0.8) +#,fill="#69b3a2", color="#e9ecef",
scale_x_continuous(breaks = seq(0,1,by=0.3))+
    ggtitle("Distribution of SNV across G4")+facet_wrap(~real_mut,scales="free")+
theme(axis.text.y=element_text(size=rel(0.8)),axis.text.x=element_text(size=rel(0.7)))) 

dev.off()




## we need to get context based study as well



test_stats<-data_total_2[, .(G4hunter_ref = mean(V9), G4hunter_alt = mean(V13),count=.N,sd=sd(V9),sd2=sd(V13)), by = .(V12,V16)]%>% fwrite(file="CLINVAR_g4hunter.tsv")


data_dfannot1%>% filter(V19=="protein-coding" &
 V8 %in%  c("promoter-TSS","exon","TTS","5' UTR","3' UTR"))%>% mutate(chr=gsub(":.*$", "", name),start=start+g4_start,end=start+g4_start+g4_stop,score=".")%>%
select(chr,start,end,name,score,G_quadstrand)%>%
fwrite(.,"G4_aroundSNVCLINVAR_all_exceptintron.bed",sep="\t",col.names=FALSE,quote=FALSE)

data_dfannot1%>% mutate(chr=gsub(":.*$", "", name),start=start+g4_start,end=start+g4_start+g4_stop,score=10)%>%
select(chr,start,end,G_quadstrand,score)%>%
fwrite(.,"G4_aroundSNVCLINVAR_inputforgreat.bed",sep="\t",col.names=FALSE,quote=FALSE)

#-------------------------------------
#this works but needs extra steps to know location of G4


filename <- "G4_aroundSNVCLINVAR_all.bed"
#system(paste0("python testMethod.py ", filename), wait = TRUE)

 system(paste0("cat ", filename ," | bedtools slop -b 2000 -i - -g ../../hg38.sizes |  sort-bed - |
 bedtools merge -s -i - -c 4,6 -o collapse,distinct -delim '+' | 
bedops --chop 100 - > G4_aroundSNVCLINVAR_flank2000_window100.bed"))

###
#sort -k1,1n -k2,2n G4_aroundSNVCLINVAR_all.bed  > G4_aroundSNVCLINVAR_all.sorted.bed
#with this we can find distance from G4 and count base

system(paste0("sortBed -i ",filename, "  | 

bedtools merge -s  -i - -c 4,6 -o collapse,distinct -delim '+'  > G4_aroundSNVCLINVAR_all.sortedmerged.bed" ))
###

system(paste0("cat G4_aroundSNVCLINVAR_flank2000_window100.bed | 
bedtools map -a - -b /bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/G4_aroundSNVCLINVAR_all.sortedmerged.bed  -c 4 -o count > G4_aroundSNVCLINVAR_all.sortedmerged.counts_G4.bed"))

system(paste0(" awk '{if ($4>0) {print} }' G4_aroundSNVCLINVAR_all.sortedmerged.counts_G4.bed > G4_region_touseCLINVAR.bed"))


##to make VCF CLINVAR FILE add chr, and sort .. need to do this only once
#cat /bio/home/goonerrn/g_quadruplex/annotation_database/SNP/ClinVar/clinvar_SNV.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "LC_ALL=C sort -k1,1 -k2,2n"}' | \
#awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' > \
# /bio/home/goonerrn/g_quadruplex/annotation_database/SNP/ClinVar/clinvar_SNV.sorted.vcf
###

system(paste0("bedtools map -a G4_aroundSNVCLINVAR_flank2000_window100.bed -b /bio/home/goonerrn/g_quadruplex/annotation_database/SNP/ClinVar/clinvar_SNV.sorted.vcf  -c 4 -o count > G4_aroundSNVCLINVAR_all.sortedmerged.counts_SNV.bed"))

##_-----we can use above two files, G4 and SNV counts to get G4 counts in each and get distance with 
#g4 file which we got above


system(paste0("bedtools closest  -D ref -t all  -a G4_aroundSNVCLINVAR_all.sortedmerged.counts_SNV.bed -b  G4_region_touseCLINVAR.bed > variants_aroundG4.clinvar.withG4closest.bed"))

#bedtools intersect  -loj -a variants_aroundG4.cosmic.bg -b G4_aroundSNVcosmic_all.sortedmerged.bed > variants_aroundG4.cosmic.withG4.bed

totalG4<-fread("G4_aroundSNVCLINVAR_all.sortedmerged.counts_G4.bed")


#flankG4SNV<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/variants_aroundG4.clinvar.withG4closest.bed")
#this is proteincoding only, and exceptintron done separately
#flankG4SNV<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/variants_aroundG4.clinvar.withG4closest_proteincoding.bed")

#flankG4SNV<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/variants_aroundG4.clinvar.withG4closestexceptintron.bed")

flankG4SNV0<-flankG4SNV%>% filter(V9==0)
flankG4SNVflanknot0<-flankG4SNV%>% filter(V9!=0)
flankG4SNVflanknot0$distgroup<-cut(flankG4SNVflanknot0$V9,seq(from=-2000,to=2000,by=100),right=FALSE,include.lowest=FALSE)
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

#temp<-split(flankSNV,flankSNV$V9)


ggplot(flankSNV, aes(x = distgroup, y = (V4),fill=distgroup)) +   geom_boxplot() +
    xlab("class") +
    theme(legend.position="none") +
    xlab("SNV in/flanking(+/-2000bp) G4")



tgc<-flankSNV%>%
  dplyr::group_by(distgroup1) %>%
  summarise(mean.count = mean(V4, na.rm = TRUE),
            sd.count = sd(V4, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se.count = sd.count / sqrt(n.count),
         lower.ci.count = mean.count - qt(1 - (0.05 / 2), n.count - 1) * se.count ,
         upper.ci.count = mean.count + qt(1 - (0.05 / 2), n.count - 1) * se.count )%>% as.data.frame()#%>% 
#as.data.frame()%>%fwrite(.,file="G4_clinvar_exceptintronpc_final.tsv",sep="\t")

pd <- position_dodge(0.1) # move them .05 to the left and right
pdf(file= "/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/pg4_clinvar_exceptintronpc_flank_count.pdf",width=10,height=7)
#names(tgc)<-c("positive strand","Negative Strand")
ggplot(tgc, aes(x=distgroup1, y=mean.count)) + 
    geom_errorbar(aes(ymin=mean.count-lower.ci.count, ymax=mean.count+upper.ci.count), colour="black", width=.2, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3)+scale_colour_Publication()+
theme_bw()+ xlab("Distance from G quadruplex divided into equal 100 bp window") + ylab("Counts of SNV")+
# theme_Publication()+
theme(axis.text.x = element_text(angle = 30, hjust = 1), 
legend.title=element_blank(), legend.position="bottom", legend.key = element_rect(color = 'white'))


dev.off()

#--------------------------------------------------------------



 a %>% mutate(sample = row_number()) %>% 
      pivot_longer(-sample) %>%
      # here you have to unselect the name column, since this is actually wrong
      select(-name) %>% 
      # separate the strings that contain the data
      separate(value, into = c("group", "value"), sep = "=") %>% 
      # put it back to wide format
      pivot_wider( names_from = group, 
                   values_from = value)




##for snpnexus
data_dfannot1%>% separate(real_mut,into=c("REF","ALT"),sep="->")%>%select(left,right,REF,ALT)%>% 
summarise(chrom=left,pos=right,ID=".",REF,ALT,QUAL="32",FILTER="PASS",INFO=".")%>% distinct()%>%
fwrite("CLINVAR_variants_inG4SNPnexusinput.vcf",sep="\t",quote=FALSE,col.names=FALSE)
##for great
data_dfannot1%>% separate(real_mut,into=c("REF","ALT"),sep="->")%>%select(left,right,REF,ALT)%>% 
summarise(chrom=left,pos=right,ID=".",REF,ALT,QUAL="32",FILTER="PASS",INFO=".")%>% distinct()%>%
fwrite("CLINVAR_variants_inG4SNPnexusinput.vcf",sep="\t",quote=FALSE,col.names=FALSE)


save(data_dfannot1,file="CLINVAR_G4_variants.RData")



##input results from snpnexus and get g4
cc_temp<-fread("cadd_g4_clinvar.txt")


cc_temp$Variant<-gsub("/","->",cc_temp$Variant)
cc_temp$Position<-as.character(cc_temp$Position)

cc_temp1<-left_join(cc_temp,data_dfannot1,by=c("Chromosome"="left","Position"="right",
"Variant"="real_mut")

cc_temp1 %>% 
  summarise(V8=V8,hours = cut(PHRED, breaks=seq(0,40, by=10), include.lowest=FALSE)) %>%
  table() %>% data.frame()# %>% arrange(desc(Freq))

cc_temp1 %>% filter(PHRED>30)%>% select(Chromosome,Position,which1,rel_loc,Variant,g4seq,G_quadstrand,V8,V19)
cc_temp1 %>% filter(PHRED>29.99)%>%
  summarise(which1=which1,hours = cut(PHRED, breaks=seq(10,40, by=10), include.lowest=FALSE)) %>%
  table() %>% data.frame()# %>% arrange(desc(Freq))



#we try to get the sequence

str_sub(data_dfannot1$G4_refrc, data_dfannot1$g4_start,  data_dfannot1$g4_stop)

str_sub(data_dfannot1$G4_refrc, data_dfannot1$g4_start,  29-data_dfannot1$g4_start)
cc_temp1%>%
mutate(a=str_sub(G4_refrc, g4_start,29),
b=str_sub(G4_refrc, 31,g4_stop))%>% head()

data_dfannot1 %>% select(G4_refrc,G4_altrc,SNP_strand,G_quadstrand)%>% head()

data_dfannot1 %>% select(G4_refrc,G4_altrc,SNP_strand,g4_start,g4_stop,G_quadstrand)%>% 
mutate(a=paste0(str_sub(G4_refrc, g4_start,29),"{",SNP_strand,"}",
str_sub(G4_refrc, 31,g4_stop)))%>% head()






#------------------------------------------------------------
load("CLINVAR_G4_variants.RData")
data_dfannot1
##get trinucleotide and pentanucleotide context

g4_maf<-data_dfannot1
library(BSgenome.Hsapiens.UCSC.hg38)

##get X nucleotide context
#g4_maf<-data_df
library(BSgenome)

genome<-BSgenome.Hsapiens.UCSC.hg38


 g4_maf$loc_end<-g4_maf$loc_start
head(g4_maf)
g4_maf<-setDT(g4_maf)
#

# we do this because the BSgen
#g4_maf<- g4_maf%>% 
#       mutate(orig_refalt = case_when((G_quadstrand=="+") ~ g4_maf$real_mut,
#                                        (G_quadstrand=="-") ~ chartr('ATGC', 'TACG',g4_maf$real_mut))) 

#for cosmic
g4_maf<-separate(g4_maf,col="real_mut",into=c("REF","ALT"),remove=FALSE)


#g4_maf$chr<-paste0("chr",g4_maf$chr)

  
g4_maf$loc_end<-g4_maf$loc_start
#g4_maf$chr<-paste0("chr",g4_maf$chr)
#g4_maf<-g4_maf[((nchar(g4_maf$ref)==1) | (nchar(g4_maf$alt)==1)),]
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
colnames(seq_right)<-"rightseq"


revString <- function(text){
  paste(rev(unlist(strsplit(text,NULL))),collapse="")
}
 

rm(seq_g4_1,seq_g4_2)
seq_g4_full<-as.data.frame(seq_g4_full)
colnames(seq_g4_full)<-"seq_g4_full"
seq<-cbind(as.data.frame(seq_left),as.data.frame(g4_maf_s),as.data.frame(seq_right),as.data.frame(seq_g4_full))

seq$G4_ref<-paste0(seq$leftseq,seq$REF,seq$rightseq)
seq$G4_alt<-paste0(seq$leftseq,seq$ALT,seq$rightseq)
                                                                             
all_1<-seq
transitions<-c("A->G","G->A","C->T","T->C")
transversions<-c("A->T","A->C","C->A","C->G","T->A","T->G","G->C","G->T")
all_1$typemut<-ifelse(all_1$real_mut %in% transitions,"Transition","Transversion")


all_2<- all_1 %>% dplyr::select(SNP_strand,real_mut,G_quadstrand,REF,ALT,which1,seq_g4_full,V8,V19,G4_ref,G4_alt,name,typemut)

      # mutate(G4_ref = case_when((G_quadstrand=="+") ~ all_1$G4_ref,
       #                                 (G_quadstrand=="-") ~ intToUtf8(rev(utf8ToInt(G4_ref)))),
	#	G4_alt=case_when((G_quadstrand=="+") ~ all_1$G4_alt,
         #                               (G_quadstrand=="-") ~ intToUtf8(rev(utf8ToInt(G4_alt)))))



temp<-split(all_2,all_2$G_quadstrand)
 map(temp,~.x%>% dplyr::group_by(G4_ref,G4_alt)%>% summarise(n=n())%>% arrange(desc(n)))%>%
rbindlist(id="source") %>% fwrite(.,file="/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/COSMIC_bychrreal_mut.tsv")


#dplyr::select(-c(leftseq,snp_left,snp_right,G4_id,rightseq))%>% rowwise()%>%
#----------------------------------------------------------------------------------------

m1<-  all_2 %>% filter(G_quadstrand=="+")%>% filter(!(is.na(G4_ref)| is.na(G4_alt)))%>%              
 dplyr::mutate(V2factor=case_when( grepl('^[G]G[ATC]GG$',G4_alt) | 
grepl('^[G][G][ATC][ATC][ATCG]$',G4_alt) | 
                                    grepl('^[ATCG][ATC][ATC][G][G]$',G4_alt) |
                                    grepl('^[ATC][G][ATC][G][ATC]$',G4_alt) | 
                                    grepl('^[ACT][G][ACT]GG$',G4_alt) | 
                                    grepl('^GG[ACT][G][ACT]$',G4_alt) ~"Nearly G4",
                                  grepl('^[ATCG][ATCG]GGG$',G4_alt) | 
grepl('^[ATCG][G][G][G][ATCG]$',G4_alt) |
                     grepl('^GG[G][ATCG][ATCG]$',G4_alt) ~"G4",       
                   TRUE ~ "loop")) %>%rowwise()%>% 
  mutate(V1factor = case_when( grepl('^[G]G[ATC]GG$',G4_ref) | grepl('^[G][G][ATC][ATC][ATCG]$',G4_ref) | 
                                  grepl('^[ATCG][ATC][ATC][G][G]$',G4_ref) |
                                  grepl('^[ATC][G][ATC][G][ATC]$',G4_ref) |
                                  grepl('^[ACT][G][ACT]GG$',G4_ref) | 
                                  grepl('^GG[ACT][G][ACT]$',G4_ref)  ~"Nearly G4",
    grepl('^[ATCG][ATCG]GGG$',G4_ref) | grepl('^[ATCG][G][G][G][ATCG]$',G4_ref) |
                           grepl('^GG[G][ATCG][ATCG]$',G4_ref) ~"G4",
                           TRUE ~ "loop")) %>% as.data.frame()



m2<-  all_2 %>% filter(G_quadstrand=="-")%>% filter(!(is.na(G4_ref)| is.na(G4_alt)))%>%              
 dplyr::mutate(V2factor=case_when( grepl('^[C]C[ATG]CC$',G4_alt) | 
grepl('^[C][C][ATG][ATG][ATGC]$',G4_alt) | 
                                    grepl('^[ATCG][ATG][ATG][C][C]$',G4_alt) |
                                    grepl('^[ATG][C][ATG][C][ATG]$',G4_alt) | 
                                    grepl('^[AGT][C][AGT]CC$',G4_alt) | 
                                    grepl('^CC[AGT][C][AGT]$',G4_alt)~"Nearly G4",
                                  grepl('^[ATCG][ATCG]CCC$',G4_alt) | 
grepl('^[ATCG][C][C][C][ATCG]$',G4_alt) |
                     grepl('^CC[C][ATCG][ATCG]$',G4_alt) ~"G4",       
                   TRUE ~ "loop")) %>%rowwise()%>% 
  mutate(V1factor = case_when( grepl('^[C]C[ATG]CC$',G4_ref) | grepl('^[C][C][ATG][ATG][ATCG]$',G4_ref) | 
                                  grepl('^[ATCG][ATG][ATG][C][C]$',G4_ref) |
                                  grepl('^[ATG][C][ATG][C][ATG]$',G4_ref) |
                                  grepl('^[AGT][C][AGT]CC$',G4_ref) | 
                                  grepl('^CC[AGT][C][AGT]$',G4_ref)  ~"Nearly G4",
    grepl('^[ATCG][ATCG]CCC$',G4_ref) | grepl('^[ATCG][C][C][C][ATCG]$',G4_ref) |
                           grepl('^CC[C][ATCG][ATCG]$',G4_ref) ~"G4",
                           TRUE ~ "loop")) %>% as.data.frame()



penclinvar<-all_2

temp<-map(temp,~setDT(.x))
map(temp,~.x[ , .N, by = .(G4_ref,G4_alt,V8) ]%>% 
dcast(paste0(G4_ref,"::",G4_alt)~V8,value.var="N")%>% arrange(desc(exon))%>% head(20)%>%as.data.frame())%>% head() #cb(.)


map(temp,~.x[ , .N, by = .(V8,which1) ]%>% 
dcast(V8~which1,value.var="N"))%>%rbindlist(id="source")#%>%cb(.)




#------------------------------------------------------------------




m1<-  all_2 %>% filter(G_quadstrand=="+")%>% filter(!(is.na(G4_ref)| is.na(G4_alt)))%>%              
  dplyr::mutate(V1factor=case_when( grepl('^[ACTG][ATC][ACTG]$',G4_ref) | 
                      grepl('^[ACGT][ACGT][ACT]$',G4_ref) |
			grepl('^[ACT][ACGT][AGCT]$',G4_ref) ~"loop",
                      grepl('^[G][G][G]$',G4_ref) ~"G4")) %>%rowwise()%>% 
  mutate(V2factor = case_when(  grepl('^G[ATC]G$',G4_alt) | 
                      grepl('^[ACT][ACGT][ACGT]$',G4_alt) |
	              grepl('^[ACGT][ACGT][ACT]$',G4_alt) ~"loop",
                      grepl('^[G][G][G]$',G4_alt) ~"G4")) %>% as.data.frame()


m2<-  all_2 %>% filter(G_quadstrand=="-")%>% filter(!(is.na(G4_ref)| is.na(G4_alt)))%>%              
 dplyr::mutate(V1factor=case_when( grepl('^[ACTG][ATG][ACTG]$',G4_ref) | 
                      grepl('^[ACGT][ACGT][AGT]$',G4_ref) |
			grepl('^[AGT][ACGT][AGCT]$',G4_ref) ~"loop",
                      grepl('^[C][C][C]$',G4_ref) ~"G4")) %>%rowwise()%>% 
  mutate(V2factor = case_when(  grepl('^C[ATG]C$',G4_alt) | 
                      grepl('^[AGT][ACGT][ACGT]$',G4_alt) |
	              grepl('^[ACGT][ACGT][AGT]$',G4_alt) ~"loop",
                      grepl('^[C][C][C]$',G4_alt) ~"G4")) %>% as.data.frame()


#-----------------------------------------------------------------------








#all_2<-all_1%>% rowwise()%>% 
#mutate(G4_ref=ifelse(all_1$G_quadstrand=="-",intToUtf8(rev(utf8ToInt(G4_ref))),
#G4_ref),
#G4_alt=ifelse(all_1$G_quadrstrand=="-",intToUtf8(rev(utf8ToInt(G4_alt))),G4_alt))

#reversed_string <- intToUtf8(rev(utf8ToInt(G4_ref)))

#triclinvar<-all_1
penclinvar<-all_2
save(triclinvar,file="/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/CONTEXT_SNVcontent_G4CLINVAR.RData")

##context study


load("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/CONTEXT_SNVcontent_G4CLINVAR.RData")
tricosmic<-all_1

head(m)
m$Var1<-as.factor(m$Var1)

#all_1<-split(tricosmic,tricosmic$G_quadstrand)
m<-table(all_1$G4_ref,all_1$G4_alt)%>% as.data.frame()%>% arrange(desc(Freq))%>% head(20)%>%
 as.data.frame()#%>% write.table(., "clipboard-16384", sep="\t", row.names=FALSE,quote=FALSE)
