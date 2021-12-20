# script to plot fisher combined transcripts fo PUS7KD and WT samples

library(tidyverse)
library(ggpubr)
library(GGally)
library(tidyr)
library(R.utils)
library(reshape2)
library(rlang)



ROOTDIR="/Volumes/scratch/TSSM/cugolini/cov"
RESULTSDIR="/Volumes/scratch/FN/TL/cugolini/cov/scripts/downstream/results"
DATADIR="/Volumes/scratch/FN/TL/cugolini/cov/analysis"
DATADIR2="/Volumes/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/nanocompore/comparison"
dir.create(RESULTSDIR)

########################### PARAMETERS ##########################################

LOR_thresh <- 1
pval_thresh <- 0.01
n_samples <- 4
share_thresh <- 4
IVT_junc_interval_left <- 25
IVT_junc_interval_right <- 25
ORF_junc_interval_left <- 15
ORF_junc_interval_right <- 15
burrows_sites <- c(22322,23317,27164,28417,28759,28927,29418)

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

ddp2 <- function(relpath){
  return(paste0(DATADIR2,"/",relpath))
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
    temp <- seq(from = i-ORF_junc_interval_left, to = i+ORF_junc_interval_right, by = 1)
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
sitelist <- read_tsv(bdp("scripts_new/backupped_data/sites_blacklist.txt"),col_types="ccc") %>% rename(`Modification type`=IUPAC) 
tx <- read_tsv(bdp("analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/orf_annotate/orf_annotate.bed"), col_names=c("chr", "start", "end", "name", "score", "strand", "cdsStart", "cdsEnd", ".", "ex", "exLen", "exSt"), col_types="cnncncnnnncc")
tx_lengths <-read.table(bdp("analysis/recappable_assembly/two_datasets/assemblies/pinfish/aln_consensus.bed"), sep = '\t',header = FALSE) %>%
  separate(V11, into=c("ex1","ex2","ex3") ,sep=",", remove=F) 
fragments <- read_tsv(bdp("scripts_new/backupped_data/RAPID/fragments_genomic_coord_UCSC.txt"),col_types = "cnn")

# CaCo2 DRS databases available
WT_C37_tx <- list.files(path = ddp("PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/WT_IVT"),pattern = "*_results.tsv" , full.names = TRUE,  recursive = T)
WT_C34_tx <- list.files(path = ddp("PUS7_KD/map_to_recap_assembly/NANOCOMPORE/sampcomp/WT"),pattern = "*_results.tsv" , full.names = TRUE,  recursive = T)
MATTHEWS_tx <- list.files(path = ddp2("matthews_caco"),pattern = "*_results.tsv" , full.names = TRUE,  recursive = T)
SRAFF_tx <- list.files(path = ddp2("sraf_caco2"),pattern = "*_results.tsv" , full.names = TRUE,  recursive = T)


##################### ADDITIONAL DATAFRAME PROCESSING ##########################

# Dataframe that returns ORFs for every transcript of the assembly
names <- select(tx, orig=name) %>% separate(orig, into=c("id", "protein"), sep="#", remove=F) %>%
  mutate(protein=case_when(is.na(protein)~"Unknown", T~protein)) %>%
  separate(protein, into=c("sp", "uniprot_id", "protein"), sep="\\|", remove=F) %>%
  select(-sp) %>%
  mutate(name=gsub("([^\\(]+).+", "\\1", protein),
         tip=as.numeric(gsub("([^\\(]+)\\(([^%]+)%/([^%]+)%\\)", "\\2", protein)),
         qip=as.numeric(gsub("([^\\(]+)\\(([^%]+)%/([^%]+)%\\)", "\\3", protein))) %>%
  select(-protein) %>%
  select(-orig) %>%
  mutate(name=case_when(is.na(name)~"Unknown", T~name))

names$name[names$id == "efad7b96-ac2e-4ce1-9b83-863ffdb18eac|116::NC_045512v2:11-29873"] <- "ORF10_SARS2"   # manually add ORF9d and ORF10 
names$name[names$id == "de81ef19-655d-4ced-a9cc-cb8384001058|107::NC_045512v2:11-29874"] <- "ORF9D_SARS2"

# Dataframe that returns transcript length for every transcript of the assembly
tx_lengths[is.na(tx_lengths)] <- 0
tx_lengths <- as.data.frame(tx_lengths) %>% mutate(sumrow= as.numeric(ex1)  + as.numeric(ex2)+as.numeric(ex3))

# Dataframe that reports IVT junctions
junction_sites <- IVT %>% 
  select(start,end) %>% 
  unlist(use.names=FALSE) %>% 
  sort()
blacklist_IVT <- vector()
for (i in junction_sites){
  temp <- seq(from = i-IVT_junc_interval_left, to = i+IVT_junc_interval_right, by = 1)
  blacklist_IVT <- c(blacklist_IVT,temp)
}
blacklist_IVT <- blacklist_IVT[which(blacklist_IVT > 0)]

# Dataframe that reports the assembly junction sites
assembly_junction_sites <- assembly %>% select(start,end,id,blockSizes,blockStarts) %>% 
  separate(blockSizes, into = c("blSize1","blSize2","blSize3"), sep = ",") %>%
  separate(blockStarts, into = c("blStart1","blStart2","blStart3"), sep = ",") %>%
  mutate(blStart1 = (as.numeric(blStart1) + as.numeric(start))) %>%
  mutate(blStart2 = (as.numeric(blStart2) + as.numeric(start))) %>%
  mutate(blStart3 = (as.numeric(blStart3) + as.numeric(start))) %>%
  mutate(blEnd1 = (as.numeric(blStart1) + as.numeric(blSize1))) %>%
  mutate(blEnd2 = (as.numeric(blStart2) + as.numeric(blSize2))) %>%
  mutate(blEnd3 = (as.numeric(blStart3) + as.numeric(blSize3))) 

# Dataframe that reports canonicity for every transcript of the assembly
canonicity <- select(assembly,start,end,id) %>%
  mutate(canonicity=ifelse(start>100,"NC", "C")) %>%
  mutate(canonicity=ifelse(end<29000,"NC", canonicity)) %>%
  mutate(canonicity=ifelse(id=="efad7b96-ac2e-4ce1-9b83-863ffdb18eac|116","NC", canonicity)) %>%  #ORF10
  mutate(canonicity=ifelse(id=="de81ef19-655d-4ced-a9cc-cb8384001058|107","NC", canonicity)) %>%  #ORF9D
  select(-start,-end) %>%
  rename(ref_id=id)
  

####################### PROCESSING OF THE DATABASES #######################


# Extraction of the datasets from Nanocompore data
WT_C37_list <- lapply(WT_C37_tx, function(x) {
  x <- read_tsv(x, col_types = "ncncccnnncncnc") %>%
    separate(cluster_counts,
             into = c("SAMPLEID", "Y", "Z"),
             sep = "_(?=[0-9])",
             remove = F
    )%>%
    select(-Y,-Z)%>%
    subset(SAMPLEID != "NC")%>%
    mutate(ref_kmer = gsub("T", "U", ref_kmer))
})


WT_C34_list <- lapply(WT_C34_tx, function(x) {
  x <- read_tsv(x, col_types = "ncncccnnncncnc") %>%
    separate(cluster_counts,
             into = c("SAMPLEID", "Y", "Z"),
             sep = "_(?=[0-9])",
             remove = F
    )%>%
    select(-Y,-Z)%>%
    subset(SAMPLEID != "NC")%>%
    mutate(ref_kmer = gsub("T", "U", ref_kmer))
})

MATTHEWS_list <- lapply(MATTHEWS_tx, function(x) {
  x <- read_tsv(x, col_types = "ncncccnnncncnc") %>%
    separate(cluster_counts,
             into = c("SAMPLEID", "Y", "Z"),
             sep = "_(?=[0-9])",
             remove = F
    )%>%
    select(-Y,-Z)%>%
    subset(SAMPLEID != "NC")%>%
    mutate(ref_kmer = gsub("T", "U", ref_kmer))
})

SRAFF_list <- lapply(SRAFF_tx, function(x) {
  x <- read_tsv(x, col_types = "ncncccnnncncnc") %>%
    separate(cluster_counts,
             into = c("SAMPLEID", "Y", "Z"),
             sep = "_(?=[0-9])",
             remove = F
    )%>%
    select(-Y,-Z)%>%
    subset(SAMPLEID != "NC")%>%
    mutate(ref_kmer = gsub("T", "U", ref_kmer))
})


# Select columns and rename samples
WT_C37 <- as.data.frame(bind_rows(WT_C37_list)) %>%
  select(-strand,-KS_dwell_pvalue,-KS_intensity_pvalue,-GMM_cov_type,-GMM_n_clust,-cluster_counts)%>%
  mutate(SAMPLEID="WT_C37")

oldnames = c("pos","genomicPos","ref_id","chr","ref_kmer","SAMPLEID")
newnames <- colnames(WT_C37)[!(colnames(WT_C37) %in% oldnames)]

WT_C37 <- WT_C37 %>%
  rename_with(~ paste0(.x, "_",unique(WT_C37$SAMPLEID)),.cols = newnames)%>%
  select(-SAMPLEID)

WT_C34_list <-WT_C34_list[sapply(WT_C34_list, function(x)dim(x)[1]) > 0]    
WT_C34 <- as.data.frame(bind_rows(WT_C34_list)) %>%
  select(-strand,-KS_dwell_pvalue,-KS_intensity_pvalue,-GMM_cov_type,-GMM_n_clust,-cluster_counts)%>%
  mutate(SAMPLEID="WT_C34")%>%
  rename_with(~ paste0(.x, "_WT_C34"),.cols = newnames)%>%
  select(-SAMPLEID)

SRAFF <- as.data.frame(bind_rows(SRAFF_list)) %>%
  select(-strand,-KS_dwell_pvalue,-KS_intensity_pvalue,-GMM_cov_type,-GMM_n_clust,-cluster_counts)%>%
  mutate(SAMPLEID="SRAFF")%>%
  rename_with(~ paste0(.x, "_SRAFF"),.cols = newnames)%>%
  select(-SAMPLEID)

MATTHEWS <- as.data.frame(bind_rows(MATTHEWS_list)) %>%
  select(-strand,-KS_dwell_pvalue,-KS_intensity_pvalue,-GMM_cov_type,-GMM_n_clust,-cluster_counts)%>%
  mutate(SAMPLEID="MATTHEWS")%>%
  rename_with(~ paste0(.x, "_MATTHEWS"),.cols = newnames)%>%
  select(-SAMPLEID) 


# Join all the samples together
oldnames<-head(oldnames,-1)
total <- full_join(WT_C37,WT_C34, by=oldnames)
total <- full_join(total,SRAFF, by=oldnames)
total <- full_join(total,MATTHEWS, by=oldnames)
total[is.na(total)] = 0
colindex <- grep("GMM_logit_pvalue", colnames(total),value=T)
samples <- sub("GMM_logit_pvalue_", "", colindex )                              # array with the name of the samples
total_split <- split(total,total$ref_id)

toplot <- lapply(X = total_split,FUN = function(x){                             # loop over the transcript models
  x <- x %>% 
    mutate(sample_presence=rowSums(select(.,starts_with("GMM_logit_pvalue"))!=0)) # count number of samples in which the site is annotated
  colindex <- grep("GMM_logit_pvalue", colnames(x),value=T)
  colindex <- sub("GMM_logit_pvalue_", "", colindex )
  for(n in colindex){                                                           # indicate for every sample using TRUE or FALSE if its LOR and pvalue are significant according to the established thresholds
    colname1<-paste0("Logit_LOR_",n)
    colname2<-paste0("GMM_logit_pvalue_",n)
    x <- x %>%
      mutate(!!n:=ifelse(!!sym(colname1)>=LOR_thresh & !!sym(colname2)<=pval_thresh,T,F))
  }
  tmp <- as.data.frame(x)%>%select(., (ncol(x)-n_samples+1):ncol(x))
  shared <- apply(tmp, MARGIN = 1, function(x) sum(x,na.rm = TRUE))             # count the samples in which the site is significant
  x<-cbind(x,shared)
  x<- x%>%
    rowwise() %>%
    mutate(ORF = get_ORF(ref_id)) %>%
    separate(ref_id, into=c("ref_id","others"), sep="::") %>%
    left_join(canonicity,by="ref_id")%>%
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
  blacklist_junctions <- junc_blacklist(unique(x$ref_id))
  x <- assign_fragment(x)
  x$fragment_ID<-factor(x$fragment_ID,levels=c("No_Fragment","Fragment1_5UTR","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10_3UTR"))
  x$sample_presence <- factor(x$sample_presence,levels=as.character(seq(from=1, to=n_samples,by=1)))
  x <- x %>% 
    mutate(IVT = ifelse(as.integer(genomicPos) %in% blacklist_junctions, "ORF junction", IVT))
  return(x)
})




# function to select max of GMM pval between the samples and take the correspondent LOR
toplot<- lapply(X = toplot,FUN = function(x){
  colindex <- grep("GMM_logit_pvalue", colnames(x))
  x<-x %>%
    rowwise() %>%
    mutate(final_pvalue=max(c_across(colindex)))
  tmp <- x %>%
    select(colindex)%>%
    rowwise() %>%
    mutate(row_max = names(.)[which.max(c_across(everything()))])%>%
    mutate(row_max=str_remove(row_max, "GMM_logit_pvalue_"))
  x<-cbind(x,tmp$row_max)%>%
    rename(row_max=`tmp$row_max`)%>%
    mutate(row_max=sub("^", "Logit_LOR_", row_max ))
  x$final_LOR<-apply(x, 1, function(y) { y[names(y)==y[names(y)==paste0("row_max")]]})
  x <- x %>% mutate(burrows_presence=ifelse(genomicPos %in% burrows_sites,T,F))
  return(x)
})


pdf(rdp("caco2_plots_per_transcript.pdf"),height=15,width=20)
lapply(toplot,function(x){
  x %>% 
    {
      ggplot(., aes(x=abs(as.numeric(final_LOR)), y=-log10(final_pvalue), color=fragment_ID, size=sample_presence)) +
        geom_point() +
        {if (nrow(subset(x,shared>=share_thresh & IVT=="No junction"))>0) ggrepel::geom_label_repel(data=filter(.,(shared>=share_thresh)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=5)}+
        scale_color_manual(breaks = c("No_Fragment","Fragment1_5UTR","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10_3UTR"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
        ggtitle(unique(x$ORF), subtitle = paste0(unique(x$ref_id)," ",unique(x$canonicity))) +
        theme_bw(22)
    }
})
dev.off()


pdf(rdp("caco2_plots_per_transcript_5p.pdf"),height=15,width=20)
lapply(toplot,function(x){
  x %>%
    subset(genomicPos<=100)%>%
    {
      ggplot(., aes(x=abs(as.numeric(final_LOR)), y=-log10(final_pvalue), size=sample_presence)) +
        geom_point() +
        {if (nrow(subset(x,shared>=share_thresh & genomicPos<=100 & IVT=="No junction"))>0) ggrepel::geom_label_repel(data=filter(., shared>=share_thresh) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=5)}+
        scale_color_manual(breaks = c("No_Fragment","Fragment1_5UTR","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10_3UTR"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
        ggtitle(unique(x$ORF), subtitle = paste0(unique(x$ref_id)," ",unique(x$canonicity))) +
        theme_bw(22)
    }
})
dev.off()


############################# WRITE TABLES #####################################

final <- lapply(toplot,function(x){
  shared_signif_mods <- nrow(subset(x,shared>=share_thresh))
  shared_and_burrows<- nrow(subset(x,shared>=share_thresh & burrows_sites==T))
  tmp<- select(x,samples)
  count_per_sample <- apply(tmp, MARGIN = 2, function(x) sum(x,na.rm=TRUE))
  table <- as.data.frame(t(count_per_sample))
  colnames(table)<-samples
  table <- table %>%
    mutate(ref_id=unique(x$ref_id),ORF=unique(x$ORF),canonicity=unique(x$canonicity))
  table <- cbind(table, shared_signif_mods,shared_and_burrows)
  return(table)
})
final <- as.data.frame(bind_rows(final))

write.table(final, rdp("caco_sites.txt"),sep="\t",quote=F,row.names=F,col.names=T)
