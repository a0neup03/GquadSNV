library(dplyr)
library(ggplot2)
library(purrr)
library(ggtext)
library(data.table)


library(dplyr)
library(ggplot2)
library(purrr)
library(ggtext)
library(data.table)
#install.packages("ggthemes")
library(ggthemes)


cb <- function(df, sep="\t", dec=",", max.size=(200*10000)){
  # Copy a data.frame to clipboard
  write.table(df, paste0("clipboard-", formatC(max.size, format="f", digits=0)), sep=sep, row.names=FALSE, dec=dec)
}



theme_Publication <- function(base_size=12, base_family="Arial") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = 14, hjust = 0.5,color="black"),
            text = element_text(family = "Arial"),
            
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(size=14, face="bold", colour = "black",angle=90,vjust =2),
            axis.title.x = element_text(size=14, face="bold", colour = "black",vjust = -0.2),
            axis.text.x = element_text(size=14, face="bold", colour = "black"), 
            legend.text = element_text(size=14,face="bold",color="black"),
            legend.title = element_text(size=14,face="bold",color="black"),
            legend.box.spacing = unit(0, "in"),
            legend.box.background = element_rect(line = 0),  
            # axis.text.y = element_text(size=12,  colour = "black"), # unbold
            axis.text.y = element_text(size=12, face="bold", colour = "black"), # bold
            strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
            strip.text.y = element_text(size = 10, face="bold", colour = "black"),
            axis.line.x = element_line(color="black", size = 0.3),
            axis.line.y = element_line(color="black", size = 0.3),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.direction = "vertical",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            panel.border=element_blank(),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ) )
}
scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33", "#117733","#44AA99", "#999933", "#888888","#88CCEE", "#CC6677", "#DDCC77", "#882255", "#332288", "#AA4499", "#661100", "#6699CC")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33", "#117733","#44AA99", "#999933", "#888888","#88CCEE", "#CC6677", "#DDCC77", "#882255", "#332288", "#AA4499", "#661100", "#6699CC")), ...)
  
}

#------------------------------------------------------------------------------------------------------------------------------------------------------------
setwd("C:/Users/a0neup03/OneDrive - University of Louisville/Documents/SNV_paper/COSMIC/")

library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)     # hg38 genome
library(BSgenome)
available.SNPs()
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
#install.packages("ggthemes")
library(ggthemes)
a<-load("C:/Users/a0neup03/OneDrive - University of Louisville/Finaldata/COSMIC_annotatrfinal.RData")
head(data_dfannot1)
data_dfannot1%>% nrow()


chain_file=  "C:/Users/a0neup03/OneDrive - University of Louisville/Documents/hg38/hg19ToHg38.over.chain.gz"


chr <- c(1:22,"X","Y","MT")

library(rtracklayer)
library(GenomicRanges)
chain_file=  "C:/Users/a0neup03/OneDrive - University of Louisville/Documents/hg38/hg19ToHg38.over.chain"
chainObject <-import.chain(chain_file)
GSM<-makeGRangesFromDataFrame(fread("C:/Users/a0neup03/Documents/G4files/GSM3003539_w15.hits.K.bed",col.names = c("chr","start","end")))

GSM_hg38 <- as.data.frame(liftOver(GSM, chainObject))

GSM_hg38%>% head()
bedA<-data_dfannot1%>% dplyr::select(left,loc_start,name,G_quadstrand)%>% 
  reframe(CHROM=left,START=loc_start,STOP=loc_start,strand=G_quadstrand,id=name)%>% distinct()%>%
  makeGRangesFromDataFrame(.,keep.extra.columns = TRUE)
GSM_hg38<-GSM_hg38%>% dplyr::select(group,seqnames,start,end,width,strand)%>%
  reframe(group=group,CHROM=seqnames,START=start,STOP=end,strand=strand,width=width)%>%  distinct()

GSM_hg38<-makeGRangesFromDataFrame(GSM_hg38,keep.extra.columns = TRUE)

# subset f1 features that overlap with f2
result <- bedA[countOverlaps(bedA, GSM_hg38, maxgap=-1) > 0]
G4_var_GSM<-as.data.frame(result)


#save(G4_var_GSM,file="COSMIC_geo_expr_g4.RData")

G4_expr<-left_join(G4_var_GSM,data_dfannot1,by=c("id"="name"),multiple="all")


G4_expr$gquad<-sapply(G4_expr$g4_type,function(x) as.integer(unlist(strsplit(x, '\\:'))[1]))

#a<-load("C:/Users/a0neup03/OneDrive - University of Louisville/Documents/G4_g4hunter_expr_results/COSMIC_geo_expr_g4_annotated.RData")

df<-G4_expr
head(df)
df$g4genestrand<-ifelse (df$annot.strand==df$strand,"same","opposite")
df$g4genestrand<-ifelse (df$annot.strand=="*","notapplicable",df$g4genestrand)
#COSMIC_annotated[as.numeric(G4_expr$right) %in% as.numeric(COSMIC_annotated$start)  ,]


nrow(df)

#df[G4_expr$id %in% df$name,]%>% select(name)%>% distinct()%>% nrow()

df$which1.ED<-ifelse(df$ED.x > df$ED.y,"ref_more","alt_more")
df$which1.ED<-ifelse(df$ED.x==df$ED.y,"no change",df$which1.ED)

transitions<-c("A->G","G->A","C->T","T->C")
transversions<-c("A->T","A->C","C->A","C->G","T->A","T->G","G->C","G->T")

df$typemut<-ifelse(df$SNP_strand %in% transitions,"Transition","Transversion")



#now we separate RNA and DNA annotations
#COSMIC_annot<-
df$annot.type<-as.character(df$annot.type)

df_a<-df %>% filter(annot.type %in% c("hg38_genes_promoters", "hg38_genes_intergenic")) %>% 
  # Keep only the rows with specific value
  dplyr::slice(rep(1:n(), each = 2)) %>% # Repeat each row twice
  mutate(annot.group = ifelse(annot.type == "hg38_genes_promoters", c("DNA","RNA"), c("DNA","RNA")))%>%as.data.frame()


df_b<- df%>% filter(!annot.type %in% c("hg38_genes_promoters", "hg38_genes_intergenic")) %>% 
  mutate(annot.group=case_when((annot.type=="hg38_lncrna_gencode")~"RNA",
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

df<-  rbind(df_b,df_a)



df<-df %>% arrange(desc(seqnames),desc(start.x),annot.gene_id) %>% group_by_at(vars(-annot.gene_id,-annot.symbol,-annot.type,-annot.p))%>% tidyr::fill(annot.gene_id,annot.symbol)%>% ungroup()

df<-setDT(df)

df.l<-split(df,df$annot.group)


df.l<-purrr::map2(df.l,names(df.l),~get_1annot(.x,.y))
(table(df.l[[2]]$annot.p1))
(table(df.l[[1]]$annot.p1))

table(df.l[[2]]$g4genestrand)

df.l[[1]]%>% arrange(desc(n))

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



#dat3<-df.l[[2]]

#name<-"RNA"
map2(df.l,names(df.l),~tests_stat(.x,.y))

library(rstatix)
library(FSA)
tests_stat<-function(dat3,name){


print(name)
if (name=="RNA") {
  dat1 <- dat3%>% filter(g4genestrand!="notapplicable")
} else {
  dat1<-dat3
}
  
  dat1$MFediff<-(dat1$MFE.y-dat1$MFE.x)
  
#res.aov<-dat1  %>% filter(annot.p!="Nearby1to5kb")%>% mutate(annot.p=as.character(annot.p))%>%
#  filter(MFediff!=0)%>%select(id,annot.p,MFediff,which1,g4genestrand,typemut,strand)%>%distinct()%>% group_by(which1,g4genestrand)%>%
#  filter(!is.na(annot.p)) %>% anova_test(  MFediff ~ annot.p,wid=id)
#print(res.aov)
#
#res.aov<-dat1       %>% filter(annot.p!="Nearby1to5kb")%>% 
#  filter(MFediff!=0)%>%select(id,annot.p,MFediff,which1,g4genestrand,typemut,strand)%>%distinct()%>% group_by(typemut,g4genestrand)%>%
#  filter(!is.na(annot.p)) %>% anova_test(  MFediff ~ annot.p,wid=id)
#res.aov


dat2<-  dat1 %>% filter(annot.p!="Nearby1to5kb")%>% 
  select(id,annot.p,MFediff,which1,MFE.x,MFE.y,g4genestrand,strand,typemut) %>% distinct() %>% filter(!is.na(annot.p))

wilcox.test(dat2$MFE.x[(dat2$which1=="Further stabilized")],dat2$MFE.y[(dat2$which1=="Further stabilized")])
t.test(dat2$MFE.x[(dat2$which1=="Further stabilized")],dat2$MFE.y[(dat2$which1=="Further stabilized")])
wilcox.test(dat2$MFE.x[(dat2$which1=="Destabilized")],dat2$MFE.y[(dat2$which1=="Destabilized")])
t.test(dat2$MFE.x[(dat2$which1=="Destabilized")],dat2$MFE.y[(dat2$which1=="Destabilized")])


#res.posthoc <- dat2 %>% filter(annot.p!="Nearby1to5kb")%>% filter(which1!="no change")%>%
#  select(id,annot.p,MFediff,which1,g4genestrand) %>% distinct() %>% filter(!is.na(annot.p)) %>% group_by(which1)%>%
#  pairwise_wilcox_test(MFediff ~ annot.p, p.adjust.method = "bonferroni", paired = F, ref.group = "PROMOTERS")

#res.posthoc 
print("1")
dat1 %>% filter(annot.p!="Nearby1to5kb")%>% 
  select(id,annot.p,MFediff,which1,g4genestrand) %>% distinct() %>% group_by(g4genestrand,which1)%>%filter(!is.na(annot.p)) %>% kruskal.test(MFediff ~ annot.p)

#dat2<-

dat2 %>% filter(annot.p!="Nearby1to5kb")%>% filter(which1!="no change")%>% 
  select(id,annot.p,MFediff,which1,g4genestrand) %>% distinct() %>% filter(!is.na(annot.p))%>% group_by(g4genestrand)%>% 
  dunn_test(MFediff~ annot.p,p.adjust.method = "bonferroni")%>% as.data.frame()%>% fwrite(.,file=paste0("COSMIC_dunntest_genestrand_annot_overall",name,".tsv"),sep='\t',col.names=T,quote=F)


dat2 %>% filter(annot.p!="Nearby1to5kb")%>%filter(which1!="no change")%>%# filter(which1=="Destabilized")%>%
  select(id,annot.p,MFediff,which1,strand,typemut) %>% distinct() %>% filter(!is.na(annot.p))%>% group_by(strand)%>%
  dunn_test(MFediff~ annot.p,p.adjust.method = "bonferroni")%>% as.data.frame()%>% fwrite(.,file=paste0("COSMIC_dunntest_strand_annot_overall",name,".tsv"),sep='\t',col.names=T,quote=F)


dat2 %>% filter(annot.p!="Nearby1to5kb")%>% filter(which1=="Destabilized")%>%
  select(id,annot.p,MFediff,which1,g4genestrand) %>% distinct() %>% filter(!is.na(annot.p))%>% group_by(g4genestrand)%>%
  dunn_test(MFediff~ annot.p,p.adjust.method = "bonferroni")%>% as.data.frame()%>% fwrite(.,file=paste0("COSMIC_dunntest_genestrand_annot_destabilized",name,".tsv"),sep='\t',col.names=T,quote=F)

dat2 %>% filter(annot.p!="Nearby1to5kb")%>% filter(which1=="Further stabilized")%>%
  select(id,annot.p,MFediff,which1,g4genestrand) %>% distinct() %>% filter(!is.na(annot.p))%>% group_by(g4genestrand)%>%
  dunn_test(MFediff~ annot.p,p.adjust.method = "bonferroni")%>% as.data.frame()%>% fwrite(.,file=paste0("COSMIC_dunntest_genestrand_annot_stabilized",name,".tsv"),sep='\t',col.names=T,quote=F)

#dunn test
dunn.test(x=abs(dat2$MFediff[(dat2$which1=="Further stabilized")]),dat2$annot.p[dat1$which1=="Further stabilized"],method = "bonferroni")
dunn.test(dat2$MFediff[(dat2$which1=="Destabilized")],dat2$annot.p[dat2$which1=="Destabilized"],method = "bonferroni")

dat2 %>% filter(annot.p!="Nearby1to5kb")%>% filter(which1=="Destabilized")%>% 
  select(id,annot.p,MFediff,which1,g4genestrand) %>% distinct() %>% filter(!is.na(annot.p))%>%group_by(g4genestrand)%>% 
  #kruskal_test(MFediff~ annot.p)
  dunn_test(MFediff~ annot.p,p.adjust.method = "bonferroni")

dat2 %>% filter(annot.p!="Nearby1to5kb")%>%# filter(which1=="Further stabilized")%>%
  group_by(g4genestrand)%>% 
  select(id,annot.p,MFediff,which1,g4genestrand) %>% distinct() %>% filter(!is.na(annot.p))%>%
  dunn_test(MFediff~ annot.p,p.adjust.method = "bonferroni")

print("2")
#dat2<-
#  dat1 %>% filter(annot.p!="Nearby1to5kb")%>% 
# select(id,annot.p,MFediff,which1,MFE.x,MFE.y,G_quadstrand) %>% distinct() %>% filter(!is.na(annot.p))

dat2%>% 
  distinct() %>% filter(!is.na(annot.p))%>% filter(which1=="Destabilized")%>%
  select(id,annot.p,MFediff,which1,g4genestrand) %>% distinct() %>% filter(!is.na(annot.p))%>% group_by(g4genestrand)%>%
  dunn_test(MFediff~ annot.p,p.adjust.method = "bonferroni")%>% as.data.frame()%>%
  fwrite(.,file=paste0("COSMIC_dunntest_genestrand_annot_destabilizd",name,".tsv"),sep='\t',col.names=T,quote=F)

dat1%>%  filter(annot.p!="Nearby1to5kb")%>% filter(which1=="Destabilized")%>%
  select(id,annot.p,MFediff,which1,strand) %>% distinct() %>%   filter(!is.na(annot.p))%>% group_by(strand)%>%
  dunn_test(MFediff~ annot.p,p.adjust.method = "bonferroni")%>% as.data.frame()%>% fwrite(.,file=paste("COSMIC_dunntest_G4strand_annot_destabilizd",name,".tsv"),sep='\t',col.names=T,quote=F)

print("3")
}
dat3<-df.l[[2]]

ggboxplot(dat3, x = "annot.p", y = "MFE.x", color = "annot.p1")
library(viridis)
# Represent it
p <- dat3 %>% mutate(MFediff=MFE.y-MFE.x)%>%filter(which1!="no change")%>%
  ggplot( aes(x=(MFediff), fill=which1)) +
  geom_histogram( color="#e9ecef", alpha=1, position = 'dodge') + scale_fill_viridis(discrete=TRUE) +
  labs(fill="")
p

table(dat1$annot.p1,dat1$genes.id)
dat3%>% mutate(MFediff=MFE.y-MFE.x)%>% mutate(MFEabs=(abs(MFediff)))%>%
  pairwise_wilcox_test(MFEabs ~ annot.type, p.adjust.method = "bonferroni")#%>% filter(p.adj<0.05)



#----------------------------------------------------------------------------------------------------------

#towards the basics



#this is the ultimate COSMIC for G4 variants without gene annotations. here the rows are repeated. 
#so to count using data_total_2 use any specific columns and non repeated rows.
#starting with how many variants affect the G4 directly or are present just upstream or downstream of the G4
x<-G4_expr%>% select(id,strand,left,V9,V12,V16,G4_refrc,G4_altrc,start.x,end.x,start.y,region,SNP_strand,MFE.x,MFE.y,CDE.x,CDE.y,CFE.x,CFE.y,G_quadstrand,affectg4,G_quadstrand)%>% distinct()
table(x$affectg4,x$G_quadstrand)
#x$affectg4<-factor(x$affectg4,levels=c("affect_g4","mut_upstream","mut_downstream"))
#x <- x[order(x$affectg4),]
x1<-x[!duplicated(x[,c('id', 'G_quadstrand')]),]
table(x1$affectg4,x1$G_quadstrand)

#from this we select the g4 affecting directly the variants.#this is what we use now.
x2<-x1%>% filter(affectg4=="affect_g4")


#not done here for cosmic and relied upon by g4 hunter, have results of quadparser to check in backup files, if needed.
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
table(x2$V12,x2$V16)


x2%>% group_by(id,left,V9,V12,V16,G4_refrc,G4_altrc,start.x,end.x,start.y,region,SNP_strand,MFE.x,MFE.y,CDE.x,CDE.y,CFE.x,CFE.y,G_quadstrand)%>%
  summarise(n=n())%>% arrange(desc(n))%>%nrow()# head()%>%as.data.frame()
#37515 variants total in COSMIC


library(tidyr)
x2<-x2%>%rowwise()%>%mutate(SNV=(strsplit(id, '[||]')[[1]][5]))%>% separate(SNV, into=c("REF","ALT"),sep="\\&")%>% as.data.frame()
#---

#======================================================================================================================================================
#======================================================================================================================================================
#======================================================================================================================================================

library(dplyr)
require(Biostrings) 
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
x2$start.y<-as.numeric(x2$start.y)
x2$chr<-as.character(x2$left)

seq_g4_1 <- BSgenome::getSeq(genome, names = x2$chr,
                             start = as.numeric(x2$start.y-1),
                             end = as.numeric(x2$start.y-1))

seq_g4_2 <-  BSgenome::getSeq(genome, names = x2$chr,
                              start = as.numeric(x2$start.y+1),
                              end =  as.numeric(x2$start.y+1))

seq_g4_full <- getSeq(genome, names = x2$chr,
                      start = as.numeric(x2$start.y-1),
                      end = as.numeric(x2$start.y+1))


head(as.data.frame(seq_g4_1))
head(as.data.frame(seq_g4_2))
head(as.data.frame(seq_g4_full))
g4_seq<-(seq_g4_full)
g4_maf <- x2 %>% mutate(G4_id = row_number())
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

#consistency not found with G4_re and mref so we will be using mref which we checked is correct and no errors. Hence problem with G4_ref and G4_alt but no time,
#and because we know mref and malt is corect we are using it.

library(stringi)
library(stringr)
#G4_refrc abd G4_altrc are both 5' to 3'. Even if G quadruplex is in the negative strand, based on the mutations we take G4 from 5' to 3'
all_1<-seq%>%mutate(mref=ifelse(strand=="+",substr(G4_refrc, 29, 31),substr(G4_refrc, 31, 33)),
                     malt=ifelse(strand=="+",substr(G4_altrc, 29, 31),substr(G4_altrc, 31, 33)))%>%
  mutate(G4_ref1=ifelse(strand=="-",stri_reverse(chartr('ATGC', 'TACG', G4_ref)),G4_ref),
         G4_alt1=ifelse(strand=="-",stri_reverse(chartr('ATGC', 'TACG', G4_alt)),G4_alt))#


#we can check if the sequences we are comparing is correct or not
all_1%>%
  filter(strand=="+")%>% select(id,strand,mref,malt,seq_g4_full,G4_ref,G4_alt,G4_refrc,G4_altrc)%>% head()

all_1%>%
  filter(strand=="-")%>% select(id,strand,mref,malt,seq_g4_full,G4_ref,G4_alt,G4_ref1,G4_alt1,G4_refrc,G4_altrc)%>%
  tail()%>% as.data.frame()


m<-table(paste0(all_1$mref,":",all_1$malt),all_1$G_quadstrand)%>% as.data.frame()%>% arrange(desc(Freq))#%>% write.table(., "clipboard-16384", sep="\t", row.names=FALSE,quote=FALSE)

m%>% head()

m<-table(as.character(paste0(all_1$mref,"->",all_1$mref)),all_1$SNP_strand,all_1$G_quadstrand)%>% as.data.frame()%>% arrange(desc(Freq))#%>% spread(Var3,Freq)%>% 
#  arrange(desc(`+`))%>% as.data.frame()%>%head(20)

m<-table(as.character(paste0(all_1$mref,"->",all_1$mref)),all_1$SNP_strand,all_1$chr)%>% as.data.frame()%>% arrange(desc(Freq))#%>% spread(Var3,Freq)%>% 

m %>% 
  group_by(Var3) %>% 
  distinct(Var1, Var2, Var3) %>%
  ungroup %>%
  summarise(pval = chisq.test(Var2, Var1)$p.value)


table(paste0((all_1$mref),"->", (all_1$malt)),all_1$SNP_strand,all_1$G_quadstrand)%>% as.data.frame()%>% arrange(desc(Freq))%>% spread(Var3,Freq)%>% 
  arrange(desc(`-`))%>% as.data.frame()%>%head(20)

table(as.character(paste0(all_1$mref,"->",all_1$malt)),all_1$SNP_strand,all_1$G_quadstrand)%>% as.data.frame()%>% arrange(desc(Freq))%>% spread(Var3,Freq)%>% 
  arrange(desc(`+`))%>% head(20)
x<-all_1%>% group_by(mref,malt,SNP_strand,G_quadstrand)%>% summarise(Freq=n())%>% 
  group_by(SNP_strand)%>% mutate(SNP_n=sum(Freq))%>% spread(G_quadstrand,Freq,fill = 0)%>% mutate(neg_s=(`-`*100)/SNP_n,pos_s=(`+`*100)/SNP_n)%>%
  arrange(desc(`-`))%>% distinct()#%>%filter(SNP_n>50)%>%head(20)



tiff("COSMIC_G4_SNV_strand.tiff", units="cm", width=8, height=5, res=300)
cairo_pdf("COSMIC_G4_SNV_strand.pdf", width=8, height=5, bg = "transparent")

#plot count of variants by strands presence in G4
all_1%>% group_by(mref,malt,SNP_strand,G_quadstrand)%>% summarise(Freq=n())%>%
  group_by(SNP_strand)%>% mutate(SNP_n=sum(Freq),gg=paste0(mref,"->",malt))%>%mutate(s=(Freq*100)/SNP_n)%>%
  ggplot(.,aes(x = SNP_strand,  y = Freq, fill = SNP_strand))+# Grouped barplot using ggplot2
  xlab("SNV")+ ylab("Frequency")+
  geom_bar(stat = "identity",
           position = "dodge")+facet_wrap(~G_quadstrand)+  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

tiff("COSMIC_G4_SNV_g4hunter.tiff", units="cm", width=8, height=5, res=300)
cairo_pdf("COSMIC_G4_SNV_g4hunter.pdf", width=8, height=5, bg = "transparent")

#plot count of variants by presence/absence of G4
all_1%>% as.data.frame()%>%select(V12,V16,id,mref,malt,SNP_strand)%>% distinct()%>%
  mutate(g4hunter=case_when((V12==-1 | V12==1) & (V16==1 | V16==-1)~"No change",
                                     ((V12==-1 | V12==1) & (V16==0))~"G4 loss", 
                                     ((V12==0) & (V16==1 | V16==-1))~"G4 Gain",
                                     TRUE~ "sth"))%>%
  group_by(g4hunter,SNP_strand)%>%  mutate(Freq=n())%>%mutate(SNP_n=sum(Freq),s=(Freq*100)/SNP_n)%>%distinct()%>%
  ggplot(.,                                      # Grouped barplot using ggplot2
         aes(x = SNP_strand,
             y = Freq,
             fill = SNP_strand)) + xlab("SNV")+ ylab("Frequency")+labs(fill='SNV') +
  geom_bar(stat = "identity",
           position = "dodge")+facet_wrap(~g4hunter)+  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

dev.off()
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

all_1$typemut<-ifelse(all_1$SNP_strand %in% transitions,"Transition","Transversion")
all_1$which1.ED<-ifelse(all_1$ED.x > all_1$ED.y,"ALT SNV stable","REF SNV stable")
all_1$which1.ED<-ifelse(all_1$ED.x==all_1$ED.y,"no change",all_1$which1.ED)
all_1$which<-ifelse(all_1$MFE.x > all_1$MFE.y,"ALT SNV stable","REF SNV stable")
all_1$which<-ifelse(all_1$MFE.x==all_1$MFE.y,"no change",all_1$which)


all_1%>% select(id,typemut,which,G_quadstrand)%>% distinct()%>%
  group_by(typemut,which,G_quadstrand)%>% mutate(Freq=n())%>% mutate(SNP_n=sum(Freq))%>%mutate(s=(Freq*100)/SNP_n)%>%
  ggplot(.,aes(x = which,  y =s, fill = typemut))+# Grouped barplot using ggplot2
  xlab("SNV")+ ylab("Frequency")+
  geom_bar(stat = "identity",
           position = "dodge")+facet_wrap(~G_quadstrand)+  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


#→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→

tiff("COSMIC_G4_3mer_bySNV.tiff", units="cm", width=16, height=10, res=300)
cairo_pdf("COSMIC_G4_3mer_bySNV.pdf", width=16, height=10, bg = "transparent")

all_1%>% select(mref,malt,SNP_strand,typemut,id)%>% distinct()%>%group_by(mref,malt,SNP_strand)%>% summarise(Freq=n(),typemut=typemut)%>% distinct()%>%arrange(desc(Freq))%>%
  group_by(SNP_strand)%>% mutate(nFreq=n(),SNP_n=sum(Freq),gg=paste0(mref,"→",malt))%>%mutate(s=(Freq*100)/SNP_n)%>%
  
ggplot(.,aes(x = reorder(as.character(gg),-Freq),  y = Freq, fill = typemut))+# Grouped barplot using ggplot2
  xlab("SNV")+ ylab("Frequency")+labs(fill='Type of SNV') +
  geom_bar(stat = "identity",
           position = "dodge")+facet_wrap(~SNP_strand,scales="free")+  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+theme(axis.text=element_text(size=12),
                                                                              axis.title=element_text(size=14,face="bold"))

dev.off()


#dont waste time on this
#this becomes if we want to check for negative stranded G4 occuring 3' direction
table(paste0(stri_reverse(all_1$mref),"→", stri_reverse(all_1$malt)),all_1$SNP_strand,all_1$G_quadstrand)%>% as.data.frame()%>% arrange(desc(Freq))%>% spread(Var3,Freq)%>% 
  arrange(desc(`-`))%>% as.data.frame()%>%head(20)
#-------------------------------------------------------------------------------------
G4_expr$which<-ifelse(G4_expr$MFE.x > G4_expr$MFE.y,"ALT stable","REF stable")
G4_expr$which<-ifelse(G4_expr$MFE.x==G4_expr$MFE.y,"no change",all_1$which)

tabletemp<- G4_expr%>% dplyr::select(SNP_strand,V9,V12,V16,V13,id,typemut,which,rel_loc)%>% distinct()%>%
  separate(SNP_strand,into=c("ref","alt"),sep="->",remove = F)

as.data.frame.matrix(table(tabletemp$ref,tabletemp$alt))%>% 
  fwrite(.,file="COSMIC_SNV_strand.tsv",sep="\t")
as.data.frame.matrix(round(prop.table(table(tabletemp$ref,tabletemp$alt))*100,2))%>% 
  fwrite(.,file="COSMIC_SNV_strand_prop.tsv",sep="\t",row.names=TRUE)
table(tabletemp$typemut,tabletemp$which)

as.data.frame.matrix(table(tabletemp$which,tabletemp$SNP_strand))%>% 
  fwrite(.,file="COSMIC_SNV_type_affectMFE_count.tsv",sep="\t")

as.data.frame.matrix(round(prop.table(table(tabletemp$typemut,tabletemp$which))*100,2))%>% 
  fwrite(.,file="COSMIC_SNV_type_affect_perc.tsv",sep="\t",row.names = TRUE) 
as.data.frame.matrix(table(tabletemp$typemut,tabletemp$which))%>% 
  fwrite(.,file="COSMIC_SNV_type_affect_count.tsv",sep="\t",row.names = TRUE) 


as.data.frame.matrix(table(paste0(tabletemp$V12),tabletemp$V16))

round(tabletemp %>% group_by(V12,V16)%>% 
        summarise(mn=n(),V9ref=mean(V9),V9sd=sd(V9),
                  V13alt=mean(V13),sdalt=sd(V13))%>% as.data.frame(),3)%>%
  fwrite(.,file="COSMIC_G4hunter_count.tsv",sep="\t") 




#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#this is slow
m<-all_1%>%filter(!(is.na(mref)| is.na(malt)))%>%    rowwise()%>%          
  dplyr::mutate(V1factor=case_when( grepl('^[ACTG][ATC][ACTG]$',mref) | 
                                      grepl('^[ACGT][ACGT][ACT]$',mref) |
                                      grepl('^[ACT][ACGT][AGCT]$',mref) ~"loop",
                                    grepl('^[G][G][G]$',mref) ~"G4")) %>% 
  mutate(V2factor = case_when(  grepl('^G[ATC]G$',malt) | 
                                  grepl('^[ACT][ACGT][ACGT]$',malt) |
                                  grepl('^[ACGT][ACGT][ACT]$',malt) ~"loop",
                                grepl('^[G][G][G]$',malt) ~"G4")) %>% as.data.frame()
#check if all ways work the same #counts say yes
library(stringr)
m <- all_1 %>% distinct()%>% 
  mutate(g4hunter=case_when((V12==-1 | V12==1) & (V16==1 | V16==-1)~"No change",
                            ((V12==-1 | V12==1) & (V16==0))~"G4 loss", 
                            ((V12==0) & (V16==1 | V16==-1))~"G4 Gain",
                            TRUE~ "sth"))%>%
  filter(!(is.na(mref) | is.na(malt))) %>%
  mutate(V1factor = ifelse(str_detect(mref, "GGG"), "G4", "loop"),
         V2factor = ifelse(str_detect(malt, "GGG"), "G4", "loop"),
         V3factor = ifelse(str_detect(mref, "(?=^[ACTG][ATC][ACTG]$)") |
                             str_detect(mref, "(?=^[ACGT][ACGT][ACT]$)") |
                             str_detect(mref, "(?=^[ACT][ACGT][AGCT]$)"), "loop", "G4"),
         V4factor = ifelse(str_detect(malt, "(?=^[ACTG][ATC][ACTG]$)") |
                             str_detect(malt, "(?=^[ACGT][ACGT][ACT]$)") |
                             str_detect(malt, "(?=^[ACT][ACGT][AGCT]$)"), "loop", "G4"))

all_1<-setDT(all_1)
all_1[, V1factor := ifelse(str_detect(G4_ref, "GGG"), "G4", "loop")]
all_1[, V2factor := ifelse(str_detect(G4_alt, "GGG"), "G4", "loop")]



as.data.frame.matrix(round(prop.table(table(m$typemut,paste0(m$V1factor,":",m$V2factor)),1)*100,2))%>%
  fwrite(.,file="COSMIC_typemut3mer_prop.tsv",sep="\t",row.names=T) 
as.data.frame.matrix((table(m$typemut,paste0(m$V1factor,":",m$V2factor))))%>%
  fwrite(.,file="COSMIC_typemut3mer_count.tsv",sep="\t",row.names=T) 

as.data.frame.matrix((table(m$which,paste0(m$V1factor,":",m$V2factor))))%>%
  fwrite(.,file="COSMIC_which3mer_count.tsv",sep="\t",row.names = T) 
as.data.frame.matrix(round(prop.table(table(m$which,paste0(m$V1factor,":",m$V2factor)),1)*100,2))%>%
  fwrite(.,file="COSMIC_which3mer_prop.tsv",sep="\t",row.names=T) 



as.data.frame.matrix((table(m$g4hunter,paste0(m$V1factor,":",m$V2factor))))%>%
  fwrite(.,file="COSMIC_g4hunter3mer_count.tsv",sep="\t",row.names = T) 
as.data.frame.matrix(round(prop.table(table(m$g4hunter,paste0(m$V1factor,":",m$V2factor)),1)*100,2))%>%
  fwrite(.,file="COSMIC_g4hunter3mer_prop.tsv",sep="\t",row.names=T) 


#-----------------------------------------------------------------------------------------------------------------------
#for this requires annotation information present in df(df.l)
#we use df.l generated below using annotatr and separate annotation for DNA and RNA
#df.l

tabletemp<-purrr::map(df.l,~.x%>%
                        dplyr::select(SNP_strand,g4genestrand,annot.p,id,typemut,which1,rel_loc)%>% distinct()%>% 
                        left_join(.,m[,c("id","mref","malt","V1factor","V2factor")],by="id"))#%>%
#separate(SNP_strand,into=c("ref","alt"),sep="->"))

purrr::map(tabletemp,~round(prop.table(table(as.character(.x$annot.p)))*100,2))
purrr::map(tabletemp,~(round(table(as.character(.x$annot.p)))))

map(tabletemp,~as.data.frame.matrix(table(as.character(.x$annot.p),.x$SNP_strand)))%>%
  bind_rows(.id="source")%>% as.data.frame()%>%
  filter(rowSums(across(where(is.numeric)))!=0)%>%
  fwrite(.,file="cosmic_annotation_SNV_count.tsv",sep="\t",row.names = TRUE)

map(tabletemp,~as.data.frame.matrix(round(prop.table(table(.x$annot.p,.x$SNP_strand),1)*100,2)))%>%
  bind_rows(.id="source")%>% as.data.frame()%>%
  filter(rowSums(across(where(is.numeric)))!=0)%>%
  fwrite(.,file="cosmic_annotation_SNV_prop.tsv",sep="\t",row.names = TRUE)

#this is again for separate annotation based on gene strand and g4 strand difference
map2(tabletemp,names(tabletemp),~onefunc_COSMIC(as.data.frame(.x),.y))


onefunc_COSMIC<-function(firstlist,DNAorRNA){
  
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
  
  fwrite(x2,file=paste0("cosmic_annotation_g4genestrand_SNV_count_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)
  fwrite(x1,file=paste0("cosmic_annotation_g4genestrand_SNV_prop_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)
  
  fwrite(y2,file=paste0("cosmic_annotation_typemut_which_count_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)
  fwrite(y1,file=paste0("cosmic_annotation_typemut_which_prop_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)
  
  fwrite(z1,file=paste0("cosmic_typemut3mer_prop_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)
  fwrite(z2,file=paste0("cosmic_typemut3mer_count_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)
  fwrite(z3,file=paste0("cosmic_which3mer_prop_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)
  fwrite(z4,file=paste0("cosmic_which3mer_count_",DNAorRNA,".tsv"),sep="\t",row.names = TRUE)
  
  return(list(y1,x1,y2,x2,z1,z2,z3,z4))
}





#======================================================================================================================================================


tiff("COSMIC_G4_SNV_effect_typmeut_percentage_overall.tiff", units="cm", width=16, height=10, res=300,compression = "lzw")
cairo_pdf("COSMIC_G4_SNV_effect_typmeut_percentage_overall.pdf", width=8, height=5, bg = "transparent")

x<-as.data.frame((table(paste0(m$g4hunter),paste0(m$V3factor,sprintf('\u2192'),m$V4factor))))

as.data.frame(prop.table(table(paste0(m$g4hunter),paste0(m$V1factor,sprintf('\u2192'),m$V2factor)),1)*100) %>% left_join(.,x,by=c("Var1","Var2"))%>%
  filter(Freq.y>0)%>% #   group_by(Var1)%>%
  arrange(desc(Freq.y)) %>%# fwrite("SNV_effect_COSMIC_trinucleotide.context.tsv")
  ggplot(aes(as.factor(Var1), Freq.x,fill=Var2)) +   labs(x="Effect of SNV on 3mer context",y="Percentage",fill="Type of SNV")+
  geom_col( position = "dodge") + 
  geom_text(aes(label=Freq.y,x=as.factor(Var1),y=Freq.x),position = position_dodge(1),vjust=-0.5)+theme_classic()+
  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+theme(legend.text =  element_text(family = "Arial"))

dev.off()

#======================================================================================================================================================


#======================================================================================================================================================



tiff("COSMIC_SNV_effect_which1_percentage.tiff",width= 16, height= 12, units="cm", res=300,compression = "lzw")
pdf("COSMIC_SNV_effect_which1_percentage.pdf", width=8, height=6)

x<-as.data.frame((table(paste0(all_1$which),paste0(all_1$V1factor,sprintf('\u2192'),all_1$V2factor))))

as.data.frame(round(prop.table(table(paste0(all_1$which),paste0(all_1$V1factor,sprintf('\u2192'),all_1$V2factor)),1)*100,2))%>% left_join(.,x,by=c("Var1","Var2"))%>%
  filter(Freq.y>0)%>% #   group_by(Var1)%>%
  arrange(desc(Freq.y)) %>%# fwrite("SNV_effect_COSMIC_trinucleotide.context.tsv")
  ggplot(aes(x=as.factor(Var2), y=Freq.x,fill=Var1)) +   labs(x="Effect of SNV on 3mer Context",y="Frequency",fill="Type of SNV")+
  geom_col( position = "dodge") + 
  geom_text(aes(label=Freq.y,x=as.factor(Var2),y=Freq.x),position = position_dodge(1),vjust=-0.5)+theme_classic()+
  scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+theme(text = element_text(family = "Times New Roman"))

dev.off()


#--------------------------------------------------------------------------------------------------------------------------------


tabletemp<-purrr::map(df.l,~.x%>% left_join(.,m[,c("id","mref","malt","V1factor","V2factor","which")],by="id")%>%
                        dplyr::select(SNP_strand,g4genestrand,annot.p,id,typemut,which,rel_loc,V1factor,V2factor,MFE.x,MFE.y,rel_loc)%>% distinct()%>% 
                        mutate(diff=(MFE.y-MFE.x),g4genestrand=recode(g4genestrand,same="Non-Template", opposite="Template"),annot.p=as.character(annot.p))%>% as.data.frame())#%>%


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
stat.test <- purrr::map(tabletemp,~.x %>% filter(g4genestrand!="notapplicable")%>%select(id,which,g4genestrand,annot.p,diff)%>%distinct()%>%
                          dplyr::filter(diff!=0)%>% mutate(annot.p=as.character(annot.p),diff=abs(diff))%>%
                   group_by(which, g4genestrand) %>% mutate(n=n())%>% 
                   pairwise_wilcox_test(formula=diff ~ annot.p,data=.) %>%
                   adjust_pvalue(method = "bonferroni") %>%
                   add_significance()%>% as.data.frame())

stat.test <- map(stat.test,~.x %>% add_xy_position(x = "annot.p", scales = "free"))

stat.test%>% rbindlist(.,idcol="source")%>%as.data.frame()%>% 
  fwrite(.,file="COSMIC_diffMFE_annotp_strandwhich_wilcox.tsv",sep="\t",quote = F)
library(ggpubr)
plots_MFE_result<-map2(tabletemp,names(tabletemp),~plot_MFE(.x,.y))

#data<-tabletemp[[2]]
#name<-names(tabletemp[[2]])
plot_MFE<-function(data,name){
  
  stat.test <- data %>% select(id,which,g4genestrand,annot.p,diff)%>%distinct()%>% filter(diff!=0 & g4genestrand!="notapplicable")%>%
    group_by(which, g4genestrand) %>% mutate(diff=abs(diff))%>% 
    rstatix::pairwise_wilcox_test(formula=diff ~ annot.p,data=.) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance()%>% as.data.frame()
  
  stat.test <- stat.test %>% add_xy_position(x = "annot.p", scales = "free")
  
  
  p<-ggplot(data = data%>% filter(diff!=0 & g4genestrand!="notapplicable"), mapping = aes(reorder(annot.p,-abs(diff)), (diff), fill=annot.p))  +
    guides(fill=F) +
    stat_summary(fun.data = quantiles_95, geom="boxplot")+ facet_grid(which~g4genestrand,drop=TRUE,scales="free")+xlab("Annotation") +
    ylab("ΔMFE ")+
    theme(strip.background = element_rect(fill = "gray90"),
          panel.border = element_rect(colour = "black", fill = NA)) +
    scale_fill_Publication()+scale_colour_Publication()+theme_Publication()  + guides(x = guide_axis(check.overlap = TRUE, angle = 45)) 
  
  p<- p +
   stat_pvalue_manual(stat.test, hide.ns = TRUE, tip.length = 0.01, inherit.aes = FALSE)+ 
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))
  
  #p<-p+ stat_compare_means(method = "anova",label = "p.format", method.args = list(test = "Tukey"))
  tiff(paste0("COSMIC_diffMFE_annotp_strandwhich_wilcox",name,".tiff"), units="cm", width=16, height=16, res=300,compression = "lzw")
  print(p)
  dev.off()
  
  cairo_pdf(paste0("COSMIC_diffMFE_annotp_strandwhich_wilcox",name,".pdf"), width=8, height=8)
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
                scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+theme(text = element_text(family = "Times New Roman")))

names(p)<-names(tabletemp)
for (i in names(tabletemp)){
  tiff(paste0("COSMIC_distribution_SNV_g4genestrand_",i,".tiff"), units="cm",width=29,height=12,res=300,compression = "lzw")
  print(p[[i]])
  dev.off()
  #pdf(paste0("COSMIC_distribution_SNV_g4genestrand",i,".pdf"), width=8,height=6)
  #print(p[[i]])
  #dev.off()
  
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
    tiff(paste0("COSMIC_SNV_g4genestrand_annot_",name,i,".tiff"), units="cm", width=24, height=12, res=300,compression = "lzw")
    #ggsave(paste0("COSMIC_SNV_g4genestrand_annotation","_",i,".svg"), width=16, height=10)
    print(result[[i]])
    dev.off()}
  
  
}





#======================================================================================================================================================
#======================================================================================================================================================
#======================================================================================================================================================
#======================================================================================================================================================

#======================================================================================================================================================
#======================================================================================================================================================
purrr::map2(df.l,names(df.l),~tritemp(.x,.y))

tritemp<-function(data,name){
  
  temp<-data%>%left_join(.,m[,c("id","mref","malt","V1factor","V2factor","which")],by="id")%>% select(id,mref,malt,g4genestrand)%>% distinct()%>%
    mutate(g4genestrand=recode(g4genestrand,same="Non-Template", opposite="Template"))%>% filter(g4genestrand!="notapplicable")%>%
    as.data.frame()#%>%
  object_count<-as.data.frame(table(paste0(temp$mref,":",temp$malt),temp$g4genestrand))%>% 
    filter(Freq>0)%>%    group_by(Var1)%>%
    arrange((Freq))
  #fwrite(.,file="tempfile.tsv",sep="\t",row.names = TRUE) 
  
  object_count$Var1<-as.character(object_count$Var1)
  object_count$Var2<-as.character(object_count$Var2)
  
  p<-object_count %>% filter(Freq>0)%>%
    mutate(group1 = if_else(Freq < 30, "Others", Var1))%>%
    group_by(group1,Var2) %>% filter(Freq>0)%>%
    summarize(avg = mean(Freq), count = n()) %>%
    ungroup() %>%
    mutate(group = if_else(group1 == "Others",
                           paste0("Others (n =", count, ")"),
                           group1)) %>%
    mutate(group = forcats::fct_reorder(group, avg)) %>% arrange(desc(count))%>%
    ggplot() + 
    geom_col(aes(group, avg,fill=Var2)) + 
    geom_text(aes(group, avg, label = round(avg, 0)),angle =90) + labs(fill="Type of Mutation")+ ylab("count")+ scale_fill_Publication()+scale_colour_Publication()+theme_Publication()+
    facet_wrap(~Var2,scales="free") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
#coord_flip()
  
  tiff(paste0("COSMIC_tri_g4genestrand_",name,".tiff"), units="cm", width=65, height=40, res=300,compression = "lzw")
  print(p)
 dev.off()
  
  cairo_pdf(paste0("COSMIC_tri_g4genestrand_",name,".pdf"), width=25, height=10)
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
plot1<-map(df.l,~.x%>%
             mutate(g4genestrand=recode(g4genestrand,same="Non-Template", opposite="Template", notapplicable="Intergenic"))%>%
            filter(!gquad %in% c(1,2,3,11,12,13))) # %>%
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
  tiff(paste0("COSMIC_G4_percentage_quartet_linegraph_g4genestrand_",names(temp)[[i]],".tiff"), units="cm", width=16, height=10, 
       res=300,compression = "lzw")
  print(temp[[i]])
  dev.off()
  cairo_pdf(paste0("COSMIC_G4_percentage_quartet_linegraph_g4genestrand_",names(temp)[[i]],".pdf"), width=25, height=10)
  print(temp[[i]])
  dev.off()
  
}





#G4_expr$gquad<-as.factor((as.integer(G4_expr$gquad)))
library(forcats)
mypalette<-c('#e6194B', '#3cb44b', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#469990', '#9A6324', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#dcbeff', '#fabed4', '#fffac8')

temp1<-map(plot1,~.x%>%dplyr::select(id,gquad,SNP_strand,g4genestrand,annot.p)%>%distinct()%>%filter(!is.na(g4genestrand) & (annot.p!="CpG Islands"))%>%
             dplyr::select(id,gquad,SNP_strand,g4genestrand)%>% distinct()%>%
             filter(gquad %in% c(4,5,6,7,8)) %>% count(SNP_strand,gquad,g4genestrand) %>%group_by(gquad,g4genestrand) %>%          # now required with changes to dplyr::count()
             mutate(prop = prop.table(n))%>%as.data.frame())

plot2<-temp1[[1]]%>%#arrange(desc(n))%>%
  ggplot(., aes( y=n, x=gquad,fill=reorder(SNP_strand,-prop),  group = reorder(SNP_strand,-prop))) + 
  geom_bar(position="stack", stat="identity")+ xlab("Number of Quartets in G quadruplex") + ylab("Counts of SNV")+ labs(color="SNV in G4")+
  facet_wrap(~g4genestrand)+theme_Publication()+
  guides(fill = guide_legend(title = "SNV in G quadruplex")) +theme(text = element_text(family = "Times New Roman"))+scale_fill_manual(name = "Type",values = mypalette)




plot2


library(ggpubr)
tiff("COSMIC_G4_percentage_quartet_linegraph_barplot_g4genestrand1.tiff", units="cm", width=25, height=15, res=300)
pdf("COSMIC_G4_percentage_quartet_linegraph_barplot_g4genestrand.pdf", width=40, height=20)

cowplot::plot_grid(temp[[1]],plot2,ncol=1,align="v")

ggarrange(temp[[1]], plot2,
          align='hv', ncol = 2,labels=c('A', 'B'),
          common.legend = T,legend="right")
dev.off()


#======================================================================================================================================================
map2(tabletemp,names(tabletemp),~plot2(.x,.y))

plot2<-function(data,name){
  tablelist<-split(data,  data$g4genestrand)
  
  
  p<-map(tablelist,~as.data.frame(round(prop.table(table(.x$annot.p,.x$which),1)*100,2)))%>% 
    bind_rows(.id="source")%>% as.data.frame()%>% 
    filter(rowSums(across(where(is.numeric)))!=0)%>% mutate(Effect=Var2)%>%mutate(Var=NULL)%>%  filter(rowSums(across(where(is.numeric)))!=0) %>% drop_na(Var1) %>%
    ggplot( aes(x = Var1, y = Freq, fill = Effect)) + facet_grid(~source, drop=TRUE, scales="free"  )+
    geom_bar(position="stack", stat="identity")+  xlab("SNV in G quadruplex region across COSMIC database") + ylab("Percentage of SNV")+
    labs(color="SNV in G4")+  theme_classic( base_family = "Arial") +
    scale_fill_Publication()+scale_colour_Publication()+theme_Publication() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  tiff(paste0("COSMIC_G4_prop_barplot_effect",name,".tiff"), units="cm", width=8, height=5, res=300)
  print(p)
  dev.off()
  pdf(paste0("COSMIC_G4_prop_barplot_effect",name,".pdf"), width=8, height=5)
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
  tiff(paste0("COSMIC_G4_count_barplot_effect",name,".tiff"), units="cm", width=8, height=5, res=300)
  print(p)
  dev.off()
  pdf(paste0("COSMIC_G4_count_barplot_effect",name,".pdf"), width=8, height=5)
  print(p)
  dev.off()
  return(p)
}




library(tidyr)
G4_expr%>% dplyr::select(seqnames,start.x,end.x,width,strand,id)%>%distinct()%>% 
  reframe(chr="Chromosome",seqnames=gsub("chr","",seqnames), start=start.x,snp=gsub("^.*\\||", "",id ))%>%separate(.,snp,"&",into=c("REF","ALT"))%>% 
  mutate("1")#%>%fwrite(.,"cosmic_snp_nexus.input.txt",col.names = F,row.names = F,quote = F,sep="\t")





#--------------------------------------------------------------------------------------------------------------


genes_cosmic_overall<-G4_expr%>%  filter(!is.na(annot.gene_id))%>% dplyr::select(annot.gene_id)%>% 
  distinct()%>% unlist(use.names=F)

genes_clinar<-G4_expr_clinvar %>%  filter(!is.na(annot.gene_id))%>% dplyr::select(annot.gene_id)%>% 
  distinct()%>% unlist(use.names=F)

cosmic_clinvar_genes<-union(genes_clinar,genes_cosmic_overall)

genes_cosmic_stabi<-G4_expr%>% filter(which1=="Further stabilized")%>%  filter(!is.na(annot.gene_id))%>% dplyr::select(annot.gene_id)%>% 
  distinct()%>% unlist(use.names=F)

genes_cosmic_g4loss<-G4_expr%>% filter((V12==1 | V12==-1) & V16==0)%>%  filter(!is.na(annot.gene_id))%>% dplyr::select(annot.gene_id)%>% 
  distinct()%>% unlist(use.names=F)
genes_cosmic_g4gain<-G4_expr%>% filter((V12==0) & (V16==1 | V16==-1))%>%  filter(!is.na(annot.gene_id))%>% dplyr::select(annot.gene_id)%>% 
  distinct()%>% unlist(use.names=F)




genes_cosmic_g4lossonly<-genes_cosmic_g4loss[!(genes_cosmic_g4loss %in% intersect(genes_cosmic_g4gain,genes_cosmic_g4loss))]
genes_cosmic_g4gainonly<-genes_cosmic_g4gain[!(genes_cosmic_g4gain %in% intersect(genes_cosmic_g4gain,genes_cosmic_g4loss))]

go_cosmic_g4gain <- gost(query = genes_cosmic_g4gainonly,
                         ordered_query = FALSE,
                         significant = TRUE, exclude_iea = TRUE, 
                         measure_underrepresentation = FALSE, evcodes = TRUE, 
                         user_threshold = 0.05, correction_method = "g_SCS", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL)


res <- as.data.frame(go_cosmic_g4gain$result)
res<-res %>% mutate(adjustedPValue=p.adjust(p_value, method = "BH", n = length(p_value)))
res <- res[order(res$source,res$adjustedPValue),]
res <- res%>% mutate(parents= vapply(parents, paste, collapse = ", ", character(1L)))
res <- res%>% mutate(evidence_codes=vapply(evidence_codes, paste, collapse = ", ", character(1L)))
res <- res%>% mutate(intersection=vapply(intersection, paste, collapse = ", ", character(1L)))

res %>%
  filter(significant=="TRUE" & source=="KEGG" & term_size>5 & adjustedPValue<0.05)%>% select(-c(intersection,evidence_codes))%>% 
  as.data.frame()%>% arrange(adjustedPValue)


