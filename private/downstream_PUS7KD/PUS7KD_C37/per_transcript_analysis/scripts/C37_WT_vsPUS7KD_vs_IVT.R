library(tidyverse)
library(ggpubr)
library(GGally)
library(tidyr)
library(R.utils)
library(reshape2)
library(dplyr)

ROOTDIR="/Volumes/scratch/TSSM/cugolini/cov"
RESULTSDIR="/Volumes/scratch/FN/TL/cugolini/cov/scripts/downstream_PUS7KD/PUS7KD_C37/per_transcript_analysis/results"
DATADIR="/Volumes/scratch/TSSM/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp"
dir.create(RESULTSDIR)
########################### FUNCTIONS ##########################################

# Function that returns full path from basedir
bdp <- function(relpath){
  return(paste0(ROOTDIR,"/",relpath))
}

# Function that returns full path from resultsdir
rdp <- function(relpath){
  return(paste0(RESULTSDIR,"/",relpath))
}

# Function that returns full path from resultsdir
ddp <- function(relpath){
  return(paste0(DATADIR,"/",relpath))
}


# Function that returns scientific notation
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

# Function that returns a vertical lab
vlab <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

# Function that returns the Fisher p-value of an array of p-values
pvalues_fisher_method = function(pvalues){
  keep = (pvalues >= 0) & (pvalues <= 1)
  pvalues[pvalues == 0] = 1e-285
  lnp = log(pvalues)
  chisq = (-2) * rowSums(lnp)
  df = 2 * length(lnp)
  fisher_pval = stats::pchisq(chisq, df, lower.tail=FALSE)
  return(fisher_pval)
}

# Function that lists the assembly junction sites for each transcript of the assembly
junc_blacklist = function(tx) {
  corr_row <- subset(assembly_junction_sites, assembly_junction_sites$id %in% tx)
  blacklist_vec <- as.numeric(as.vector(corr_row[1,7:12]))
  blacklist_vec <- blacklist_vec[!is.na(blacklist_vec)]
  black_int <- vector()
  for (i in blacklist_vec){
    temp <- seq(from = i-15, to = i+15, by = 1)
    black_int <- c(black_int,temp)
  }
  black_int <- black_int[which(black_int > 0)] %>% sort()
  return(black_int)
}           

# Function to get a transcript length
get_tx_length = function(ref_id) {
  tx_length = tx_lengths[with(tx_lengths, tx_lengths$V4 %in% ref_id),]$sumrow
  return(tx_length)
}

# Function to get the ORF encoded in a transcript
get_ORF = function(ref_id) {
  temp <- names[with(names, id %in% ref_id),]$name
  return(temp)
}

# Function to assign fragment id to each position
assign_fragment = function(df) {
  df$fragment_ID <- "No_Fragment"
  df$genomicPos <- as.numeric(df$genomicPos)
  for (i in 1:length(fragments$ID)){
    df$fragment_ID <- ifelse(df$genomicPos >= fragments$start[i] & df$genomicPos <= fragments$end[i], fragments$ID[i], df$fragment_ID)
  }
  return(df)
}




########################### GENERAL DATA #######################################

IVT <-read_tsv(bdp("scripts_new/backupped_data/IVT_junctions.bed"), col_names = c("chrom","start","end","id","score","strand"),col_types="cnncnc")
assembly <-read.table(bdp("scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus.bed"), col.names = c("chrom","start","end","id","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"), sep="\t")
sitelist <- read_tsv(bdp("scripts_new/backupped_data/sites_blacklist.txt"),col_types="ccc") %>% dplyr::rename(`Modification type`=IUPAC) 
tx <- read_tsv(bdp("analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/orf_annotate/orf_annotate.bed"), col_names=c("chr", "start", "end", "name", "score", "strand", "cdsStart", "cdsEnd", ".", "ex", "exLen", "exSt"), col_types="cnncncnnnncc")
tx_lengths <-read.table(bdp("analysis/recappable_assembly/two_datasets/assemblies/pinfish/aln_consensus.bed"), sep = '\t',header = FALSE) %>%
  separate(V11, into=c("ex1","ex2","ex3") ,sep=",", remove=F) 
fragments <- read_tsv(bdp("scripts_new/backupped_data/RAPID/fragments_genomic_coord_UCSC.txt"),col_types = "cnn")

PUS7_KD_WT_tx <- list.files(path = ddp("PUS7_KD_WT/"),pattern = "*_results.tsv" , full.names = TRUE,  recursive = T)
PUS7_KD_IVT_tx <- list.files(path = ddp("PUS7_KD_IVT"),pattern = "*_results.tsv" , full.names = TRUE,  recursive = T)
WT_IVT_tx <- list.files(path = ddp("WT_IVT"),pattern = "*_results.tsv" , full.names = TRUE,  recursive = T)


##################### ADDITIONAL DATAFRAME PROCESSING ##########################

# Dataframe that returns ORFs for every transcript of the assembly
names <- dplyr::select(tx, orig=name) %>% separate(orig, into=c("id", "protein"), sep="#", remove=F) %>%
  mutate(protein=case_when(is.na(protein)~"Unknown", T~protein)) %>%
  separate(protein, into=c("sp", "uniprot_id", "protein"), sep="\\|", remove=F) %>%
  dplyr::select(-sp) %>%
  mutate(name=gsub("([^\\(]+).+", "\\1", protein),
         tip=as.numeric(gsub("([^\\(]+)\\(([^%]+)%/([^%]+)%\\)", "\\2", protein)),
         qip=as.numeric(gsub("([^\\(]+)\\(([^%]+)%/([^%]+)%\\)", "\\3", protein))) %>%
  dplyr::select(-protein) %>%
  dplyr::select(-orig) %>%
  mutate(name=case_when(is.na(name)~"Unknown", T~name))

names$name[names$id == "efad7b96-ac2e-4ce1-9b83-863ffdb18eac|116::NC_045512v2:11-29873"] <- "ORF10_SARS2"   # manually add ORF9d and ORF10 
names$name[names$id == "de81ef19-655d-4ced-a9cc-cb8384001058|107::NC_045512v2:11-29874"] <- "ORF9D_SARS2"

# Dataframe that returns transcript length for every transcript of the assembly
tx_lengths[is.na(tx_lengths)] <- 0
tx_lengths <- as.data.frame(tx_lengths) %>% mutate(sumrow= as.numeric(ex1)  + as.numeric(ex2)+as.numeric(ex3))

# Dataframe that reports IVT junctions
junction_sites <- IVT %>% 
  dplyr::select(start,end) %>% 
  unlist(use.names=FALSE) %>% 
  sort()
blacklist_IVT <- vector()
for (i in junction_sites){
  temp <- seq(from = i-25, to = i+25, by = 1)
  blacklist_IVT <- c(blacklist_IVT,temp)
}
blacklist_IVT <- blacklist_IVT[which(blacklist_IVT > 0)]

# Dataframe that reports the assembly junction sites
assembly_junction_sites <- assembly %>% dplyr::select(start,end,id,blockSizes,blockStarts) %>% 
  separate(blockSizes, into = c("blSize1","blSize2","blSize3"), sep = ",") %>%
  separate(blockStarts, into = c("blStart1","blStart2","blStart3"), sep = ",") %>%
  mutate(blStart1 = (as.numeric(blStart1) + as.numeric(start))) %>%
  mutate(blStart2 = (as.numeric(blStart2) + as.numeric(start))) %>%
  mutate(blStart3 = (as.numeric(blStart3) + as.numeric(start))) %>%
  mutate(blEnd1 = (as.numeric(blStart1) + as.numeric(blSize1))) %>%
  mutate(blEnd2 = (as.numeric(blStart2) + as.numeric(blSize2))) %>%
  mutate(blEnd3 = (as.numeric(blStart3) + as.numeric(blSize3))) 

# Dataframe that reports canonicity for every transcript of the assembly
canonicity <- dplyr::select(assembly,start,end,id) %>%
  mutate(canonicity=ifelse(start>100,"NC", "C")) %>%
  mutate(canonicity=ifelse(end<29000,"NC", canonicity)) %>%
  dplyr::select(-start,-end) %>%
  dplyr::rename(ref_id=id)



############################ PROCESS ###########################################



##### PUS7_KD VS WT modifications ####


PUS7_KD_WT <- lapply(PUS7_KD_WT_tx, function(x) {
  x <- read_tsv(x, col_types = "ncncccnnncncnc") %>%
    separate(cluster_counts,
      into = c("SAMPLEID", "Y", "Z"),
      sep = "_(?=[0-9])",
      remove = F
    )%>%
    dplyr::select(-Y,-Z)%>%
    subset(SAMPLEID != "NC")%>%
    mutate(ref_kmer = gsub("T", "U", ref_kmer))
})

PUS7_KD_WT <-
  PUS7_KD_WT[sapply(PUS7_KD_WT, function(x)
    dim(x)[1]) > 0]                                                             #delete 0-length dataframes
total <- bind_rows(PUS7_KD_WT) %>%
  rowwise() %>%
  mutate(ORF = get_ORF(ref_id)) %>%
  separate(ref_id, into=c("ref_id","others"), sep="::") %>%
  rowwise() %>%
  mutate(tx_length = get_tx_length(ref_id)) %>%
  mutate(start_end_sites = ifelse(pos < 100, "pstart", NA)) %>%
  mutate(start_end_sites = ifelse(pos > (tx_length - 40), "pend", start_end_sites)) %>%
  left_join(sitelist, by = "ref_kmer") %>%
  mutate(IVT = ifelse(
    genomicPos %in% blacklist_IVT,
    "IVT junction",
    "No junction"
  )) %>%
  mutate(`Modification type`=ifelse(is.na(`Modification type`),"Others",`Modification type`))
total_split <- split(total, total$ref_id)
total_split <- lapply(total_split,function(x){
  blacklist_junctions <- junc_blacklist(unique(x$ref_id))
  x <- assign_fragment(x)
  x$fragment_ID<-factor(x$fragment_ID,levels=c("No_Fragment","Fragment1_5UTR","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10_3UTR"))
  x <- x %>% 
    mutate(IVT = ifelse(as.integer(genomicPos) %in% blacklist_junctions, "ORF junction", IVT))
  return(x)
})


pdf(rdp("WT_vs_PUS7_KD_per transcript.pdf"),height=15,width=20)
lapply(total_split,function(x){
  x %>% 
    {
      ggplot(., aes(x=abs(Logit_LOR), y=-log10(GMM_logit_pvalue),  color=fragment_ID)) +
        geom_point(size=2) +
        {if (nrow(subset(x,-log10(GMM_logit_pvalue)>=1 & abs(Logit_LOR)>=0.5))>0) ggrepel::geom_label_repel(data=filter(., (-log10(GMM_logit_pvalue)>=1 & abs(Logit_LOR)>=0.5)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=5)}+
        scale_color_manual(breaks = c("No_Fragment","Fragment1_5UTR","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10_3UTR"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
        ggtitle(unique(x$ORF), subtitle = unique(x$ref_id)) +
        theme_bw(22)
    }
})
dev.off()

##  site selection
PUS7_KD_WT_selected_sites <- bind_rows(total_split) %>%
  arrange(GMM_logit_pvalue)
write.table(WT_IVT_selected_sites, rdp("PUS7_KD_WT_selected_sites.txt"), row.names = F,quote = F,sep="\t")

# I select the top 4 modifications and write correspondance tables 
CCAAU_24464 <- PUS7_KD_WT_selected_sites %>% unite(ref_id,ref_id,others, sep = "::")%>% subset(genomicPos==24464)%>% select(pos,ref_id)
write.table(CCAAU_24464,sep = "\t",quote = F,row.names = F,col.names = F,file =ddp("corresp_tx_genome_tables/CCAAU_24464.txt"))

UAAUA_23510 <- PUS7_KD_WT_selected_sites%>% subset(genomicPos==23510) %>% select(pos,ref_id)
write.table(UAAUA_23510,sep = "\t",quote = F,row.names = F,col.names = F,file =ddp("corresp_tx_genome_tables/UAAUA_23510.txt"))

GCCAU_26337 <- PUS7_KD_WT_selected_sites %>% subset(genomicPos==26337) %>% select(pos,ref_id)
write.table(GCCAU_26337,sep = "\t",quote = F,row.names = F,col.names = F,file =ddp("corresp_tx_genome_tables/GCCAU_26337.txt"))

CUUGA_28960 <- PUS7_KD_WT_selected_sites %>% subset(genomicPos==28960) %>% select(pos,ref_id)
write.table(CUUGA_28960,sep = "\t",quote = F,row.names = F,col.names = F,file =ddp("corresp_tx_genome_tables/CUUGA_28960.txt"))



##  plot only sites with genomic Position<=100

pdf(rdp("WT_vs_PUS7_KD_per transcript_5p.pdf"),height=15,width=20)
lapply(total_split,function(x){
  x<-subset(x,as.integer(genomicPos)<=100)
  if(nrow(x)>0){
    x %>% 
      {
        ggplot(., aes(x=abs(Logit_LOR), y=-log10(GMM_logit_pvalue),  color=`Modification type`)) +
          geom_point(size=2) +
          scale_color_manual(values=c("#FF0000","#999999", "#56B4E9")) +
          {if (nrow(subset(x,-log10(GMM_logit_pvalue)>=1 & abs(Logit_LOR)>=0.5))>0) ggrepel::geom_label_repel(data=filter(., (-log10(GMM_logit_pvalue)>=1 & abs(Logit_LOR)>=0.5)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=5)}+
          ggtitle(unique(x$ORF), subtitle = unique(x$ref_id)) +
          theme_bw(22)
      }
  }
})
dev.off()


##### PUS7_KD VS IVT modifications ####

PUS7_KD_IVT <- lapply(PUS7_KD_IVT_tx, function(x) {
  x <- read_tsv(x, col_types = "ncncccnnncncnc") %>%
    separate(cluster_counts,
             into = c("SAMPLEID", "Y", "Z"),
             sep = "_(?=[0-9])",
             remove = F
    )%>%
    dplyr::select(-Y,-Z)%>%
    subset(SAMPLEID != "NC")%>%
    mutate(ref_kmer = gsub("T", "U", ref_kmer))
})

PUS7_KD_IVT <-
  PUS7_KD_IVT[sapply(PUS7_KD_IVT, function(x)
    dim(x)[1]) > 0]                                                             #delete 0-length dataframes
total <- bind_rows(PUS7_KD_IVT) %>%
  rowwise() %>%
  mutate(ORF = get_ORF(ref_id)) %>%
  separate(ref_id, into=c("ref_id","others"), sep="::") %>%
  rowwise() %>%
  mutate(tx_length = get_tx_length(ref_id)) %>%
  mutate(start_end_sites = ifelse(pos < 100, "pstart", NA)) %>%
  mutate(start_end_sites = ifelse(pos > (tx_length - 40), "pend", start_end_sites)) %>%
  left_join(sitelist, by = "ref_kmer") %>%
  mutate(IVT = ifelse(
    genomicPos %in% blacklist_IVT,
    "IVT junction",
    "No junction"
  )) %>%
  mutate(`Modification type`=ifelse(is.na(`Modification type`),"Others",`Modification type`))
total_split <- split(total, total$ref_id)
total_split <- lapply(total_split,function(x){
  blacklist_junctions <- junc_blacklist(unique(x$ref_id))
  x <- assign_fragment(x)
  x$fragment_ID<-factor(x$fragment_ID,levels=c("No_Fragment","Fragment1_5UTR","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10_3UTR"))
  x <- x %>% 
    mutate(IVT = ifelse(as.integer(genomicPos) %in% blacklist_junctions, "ORF junction", IVT))
  return(x)
})



pdf(rdp("PUS7_KD_vs_IVT_per transcript.pdf"),height=15,width=20)
lapply(total_split,function(x){
  x %>% 
    {
      ggplot(., aes(x=abs(Logit_LOR), y=-log10(GMM_logit_pvalue),  color=fragment_ID)) +
        geom_point(size=2) +
        {if (nrow(subset(x,-log10(GMM_logit_pvalue)>=10 & abs(Logit_LOR)>=1))>0) ggrepel::geom_label_repel(data=filter(., (-log10(GMM_logit_pvalue)>=10 & abs(Logit_LOR)>=1)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=5)}+
        scale_color_manual(breaks = c("No_Fragment","Fragment1_5UTR","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10_3UTR"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
        ggtitle(unique(x$ORF), subtitle = unique(x$ref_id)) +
        theme_bw(22)
    }
})
dev.off()

##  plot only sites with genomic Position<=100

pdf(rdp("PUS7_KD_vs_IVT_per transcript_5p.pdf"),height=15,width=20)
lapply(total_split,function(x){
  x<-subset(x,as.integer(genomicPos)<=100)
  if(nrow(x)>0){
    x %>% 
      {
        ggplot(., aes(x=abs(Logit_LOR), y=-log10(GMM_logit_pvalue),  color=`Modification type`)) +
          geom_point(size=2) +
          scale_color_manual(values=c("#FF0000","#999999", "#56B4E9")) +
          {if (nrow(subset(x,-log10(GMM_logit_pvalue)>=1 & abs(Logit_LOR)>=1))>0) ggrepel::geom_label_repel(data=filter(., (-log10(GMM_logit_pvalue)>=1 & abs(Logit_LOR)>=1)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=5)}+
          ggtitle(unique(x$ORF), subtitle = unique(x$ref_id)) +
          theme_bw(20)
      }
  }
})
dev.off()

##  site selection

PUS7_KD_IVT_selected_sites <- bind_rows(total_split) %>%
  subset(abs(Logit_LOR)>=1 & GMM_logit_pvalue<=0.05) %>%
  arrange(GMM_logit_pvalue)%>%
  subset(fragment_ID!="No_Fragment")%>%
  mutate(third=substring(ref_kmer, 3, 3))%>%
  subset(third=="U")%>%
  left_join(canonicity,by="ref_id") %>%
  subset(canonicity=="C")
write.table(PUS7_KD_IVT_selected_sites, rdp("PUS7_KD_IVT_selected_sites.txt"), row.names = F,quote = F,sep="\t")



##### WT VS IVT modifications ####

WT_IVT <- lapply(WT_IVT_tx, function(x) {
  x <- read_tsv(x, col_types = "ncncccnnncncnc") %>%
    separate(cluster_counts,
             into = c("SAMPLEID", "Y", "Z"),
             sep = "_(?=[0-9])",
             remove = F
    )%>%
    dplyr::select(-Y,-Z)%>%
    subset(SAMPLEID != "NC")%>%
    mutate(ref_kmer = gsub("T", "U", ref_kmer))
})

WT_IVT <-
  WT_IVT[sapply(WT_IVT, function(x)
    dim(x)[1]) > 0]                                                             #delete 0-length dataframes
total <- bind_rows(WT_IVT) %>%
  rowwise() %>%
  mutate(ORF = get_ORF(ref_id)) %>%
  separate(ref_id, into=c("ref_id","others"), sep="::") %>%
  rowwise() %>%
  mutate(tx_length = get_tx_length(ref_id)) %>%
  mutate(start_end_sites = ifelse(pos < 100, "pstart", NA)) %>%
  mutate(start_end_sites = ifelse(pos > (tx_length - 40), "pend", start_end_sites)) %>%
  left_join(sitelist, by = "ref_kmer") %>%
  mutate(IVT = ifelse(
    genomicPos %in% blacklist_IVT,
    "IVT junction",
    "No junction"
  )) %>%
  mutate(`Modification type`=ifelse(is.na(`Modification type`),"Others",`Modification type`))
total_split <- split(total, total$ref_id)
total_split <- lapply(total_split,function(x){
  blacklist_junctions <- junc_blacklist(unique(x$ref_id))
  x <- assign_fragment(x)
  x$fragment_ID<-factor(x$fragment_ID,levels=c("No_Fragment","Fragment1_5UTR","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10_3UTR"))
  x <- x %>% 
    mutate(IVT = ifelse(as.integer(genomicPos) %in% blacklist_junctions, "ORF junction", IVT))
  return(x)
})



pdf(rdp("WT_vs_IVT_per transcript.pdf"),height=15,width=20)
lapply(total_split,function(x){
  x %>% 
    {
      ggplot(., aes(x=abs(Logit_LOR), y=-log10(GMM_logit_pvalue),  color=fragment_ID)) +
        geom_point(size=2) +
        {if (nrow(subset(x,-log10(GMM_logit_pvalue)>=10 & abs(Logit_LOR)>=1))>0) ggrepel::geom_label_repel(data=filter(., (-log10(GMM_logit_pvalue)>=10 & abs(Logit_LOR)>=1)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=5)}+
        scale_color_manual(breaks = c("No_Fragment","Fragment1_5UTR","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10_3UTR"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
        ggtitle(unique(x$ORF), subtitle = unique(x$ref_id)) +
        theme_bw(22)
    }
})
dev.off()

##  plot only sites with genomic Position<=100

pdf(rdp("WT_vs_IVT_per transcript_5p.pdf"),height=15,width=20)
lapply(total_split,function(x){
  x<-subset(x,as.integer(genomicPos)<=100)
  if(nrow(x)>0){
    x %>% 
      {
        ggplot(., aes(x=abs(Logit_LOR), y=-log10(GMM_logit_pvalue),  color=`Modification type`)) +
          geom_point(size=2) +
          scale_color_manual(values=c("#FF0000","#999999", "#56B4E9")) +
          {if (nrow(subset(x,-log10(GMM_logit_pvalue)>=1 & abs(Logit_LOR)>=1))>0) ggrepel::geom_label_repel(data=filter(., (-log10(GMM_logit_pvalue)>=1 & abs(Logit_LOR)>=1)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=5)}+
          ggtitle(unique(x$ORF), subtitle = unique(x$ref_id)) +
          theme_bw(20)
      }
  }
})
dev.off()

##  site selection

WT_IVT_selected_sites <- bind_rows(total_split) %>%
  subset(abs(Logit_LOR)>=1 & GMM_logit_pvalue<=0.05) %>%
  arrange(GMM_logit_pvalue)%>%
  subset(fragment_ID!="No_Fragment")%>%
  mutate(third=substring(ref_kmer, 3, 3))%>%
  subset(third=="U")%>%
  left_join(canonicity,by="ref_id") %>%
  subset(canonicity=="C")
write.table(WT_IVT_selected_sites, rdp("WT_IVT_selected_sites.txt"), row.names = F,quote = F,sep="\t")
