library(tidyverse)
library(dplyr)
library(ggpubr)
library(GGally)
library(tidyr)
library(R.utils)
library(reshape2)
library(rlang)
library(xlsx)



########################### PARAMETERS ##########################################

LOR_thresh <- 0.5
pval_thresh <- 0.01
n_samples <- 1
IVT_junc_interval_left <- 25
IVT_junc_interval_right <- 25
ORF_junc_interval_left <- 15
ORF_junc_interval_right <- 15
fleming_paper_notation <- c(22322,23317,27164,28417,28759,28927,29418)
fleming_sites_bed <- fleming_paper_notation-3
which_nucl <- "U"                                                               # possibilities are U or A


########################### DIRECTORIES ##########################################

ROOTDIR="/Volumes/scratch/TSSM/cugolini/cov"
RESULTSDIR="/Volumes/scratch/FN/TL/cugolini/cov/scripts/downstream_PUS7KD/PUS7KD_C37/per_nucleotide_analysis"
DATADIR="/Volumes/scratch/TSSM/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp"
dir.create(RESULTSDIR,showWarnings = F)


########################### FUNCTIONS ##########################################

# Function that returns full path from basedir
bdp <- function(relpath){
  return(paste0(ROOTDIR,"/",relpath))
}

# Function that returns full path from resultsdir
rdp <- function(relpath){
  return(paste0(RESULTSDIR,"/",relpath))
}

# Function that returns full path from datadir
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


# Function that lists the assembly junction sites for each transcript of the assembly
junc_blacklist = function(transcript) {
  corr_row <- subset(assembly_junction_sites, assembly_junction_sites$id %in% transcript)
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

# Function that blacklists the whole kmer for ORF and IVT junctions
blacklist_kmer = function(pos_first_nucl,transcript){
  blacklist_ORF <-junc_blacklist(transcript)
  seq_kmer <- seq(pos_first_nucl,(pos_first_nucl+4),1)
  junction_type <- ifelse(any(seq_kmer %in% blacklist_IVT)==T,"IVT junction","No junction")
  junction_type <- ifelse(any(seq_kmer %in% blacklist_ORF)==T,"ORF junction",junction_type)
  return(junction_type)
}

# Function that blacklists the whole kmer for start and end sites of transcripts
blacklist_start_end = function(tx_pos,length_tx){
  start_end_type <- ifelse(tx_pos < 100, "pstart", NA)
  start_end_type <- ifelse(tx_pos > (length_tx-40), "pend", start_end_type)
  return(start_end_type)
}

# Function to know if fleming sites are present
get_fleming = function(site) {
  blacklist_fleming <- vector()
  for (i in fleming_sites_bed){
    temp <- seq(from = i, to = i+4, by = 1)
    blacklist_fleming <- c(blacklist_fleming,temp)
  }
  return(site %in% blacklist_fleming)
}

# Function to get a transcript length
get_tx_length = function(transcript) {
  tx_length = tx_lengths[with(tx_lengths, tx_lengths$V4 %in% transcript),]$sumrow
  return(tx_length)
}

# Function to get the ORF encoded in a transcript
get_ORF = function(transcript) {
  temp <- names[with(names, ref_id %in% transcript),]$ORF
  return(temp)
}

# Function to assign fragment id to each position
assign_fragment_kmer_first_nucl = function(pos_first_nucl) {
  fragment_ID <- "No_Fragment"
  for (i in 1:length(fragments_bed$ID)){
    fragment_ID <- ifelse(pos_first_nucl >= fragments_bed$start[i] & (pos_first_nucl+4) <= fragments_bed$end[i], fragments_bed$ID[i], fragment_ID)
  }
  return(fragment_ID)
}



########################### GENERAL DATA #######################################

IVT <- read_tsv(bdp("scripts_new/backupped_data/IVT_junctions.bed"), col_names = c("chrom","start","end","id","score","strand"),col_types="cnncnc") %>% 
  dplyr::mutate(start=(start-1),end=(end-1))
assembly <- read.table(bdp("scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus.bed"), col.names = c("chrom","start","end","id","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"), sep="\t")
sitelist <- read_tsv(bdp("scripts_new/backupped_data/sites_blacklist.txt"),col_types="ccc") %>%
  dplyr::rename(motif=IUPAC) 
tx <- read_tsv(bdp("analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/orf_annotate/orf_annotate.bed"), col_names=c("chr", "start", "end", "ref_id", "score", "strand", "cdsStart", "cdsEnd", ".", "ex", "exLen", "exSt"), col_types="cnncncnnnncc")
tx_lengths <- read.table(bdp("analysis/recappable_assembly/two_datasets/assemblies/pinfish/aln_consensus.bed"), sep = '\t',header = FALSE) %>%
  separate(V11, into=c("ex1","ex2","ex3") ,sep=",", remove=F) 
fragments <- read_tsv(bdp("scripts_new/backupped_data/RAPID/fragments_genomic_coord_UCSC.txt"),col_types = "cnn")
fragments_bed <- fragments %>% 
  dplyr::mutate(start=(start-1),end=(end-1))

# DRS databases available
PUS7_KD_WT_tx <- list.files(path = ddp("PUS7_KD_WT/"),pattern = "*_results.tsv" , full.names = TRUE,  recursive = T)
WT_IVT_tx <- list.files(path = ddp("WT_IVT/"),pattern = "*_results.tsv" , full.names = TRUE,  recursive = T)



##################### ADDITIONAL DATAFRAME PROCESSING ##########################

# Dataframe that returns ORFs for every transcript of the assembly
names <- dplyr::select(tx, orig=ref_id) %>% 
  separate(orig, into=c("ref_id", "protein"), sep="#", remove=F) %>%
  mutate(protein=case_when(is.na(protein)~"Unknown", T~protein)) %>%
  separate(protein, into=c("sp", "uniprot_id", "protein"), sep="\\|", remove=F) %>%
  dplyr::select(-sp) %>%
  mutate(ORF=gsub("([^\\(]+).+", "\\1", protein),
         tip=as.numeric(gsub("([^\\(]+)\\(([^%]+)%/([^%]+)%\\)", "\\2", protein)),
         qip=as.numeric(gsub("([^\\(]+)\\(([^%]+)%/([^%]+)%\\)", "\\3", protein))) %>%
  dplyr::select(-protein) %>%
  dplyr::select(-orig) %>%
  mutate(ORF=case_when(is.na(ORF)~"Unknown", T~ORF))

names$ORF[names$ref_id == "efad7b96-ac2e-4ce1-9b83-863ffdb18eac|116::NC_045512v2:11-29873"] <- "ORF10_SARS2"   # manually add ORF9d and ORF10 
names$ORF[names$ref_id == "de81ef19-655d-4ced-a9cc-cb8384001058|107::NC_045512v2:11-29874"] <- "ORF9D_SARS2"

# Dataframe that returns transcript length for every transcript of the assembly
tx_lengths[is.na(tx_lengths)] <- 0
tx_lengths <- as.data.frame(tx_lengths) %>% 
  mutate(sumrow= as.numeric(ex1) + as.numeric(ex2)+as.numeric(ex3))

# Dataframe that reports IVT junctions
junction_sites <- IVT %>% 
  dplyr::select(start,end) %>% 
  unlist(use.names=FALSE) %>% 
  sort()
blacklist_IVT <- vector()
for (i in junction_sites){
  temp <- seq(from = i-IVT_junc_interval_left, to = i+IVT_junc_interval_right, by = 1)
  blacklist_IVT <- c(blacklist_IVT,temp)
}
blacklist_IVT <- blacklist_IVT[which(blacklist_IVT >= 0)]

# Dataframe that reports the assembly junction sites
assembly_junction_sites <- assembly %>% 
  dplyr::select(start,end,id,blockSizes,blockStarts) %>% 
  separate(blockSizes, into = c("blSize1","blSize2","blSize3"), sep = ",",convert=T) %>%
  separate(blockStarts, into = c("blStart1","blStart2","blStart3"), sep = ",",convert=T) %>%
  mutate(blStart1 = (blStart1 + start)) %>%
  mutate(blStart2 = (blStart2 + start)) %>%
  mutate(blStart3 = (blStart3 + start)) %>%
  mutate(blEnd1 = (blStart1 + blSize1)) %>%
  mutate(blEnd2 = (blStart2 + blSize2)) %>%
  mutate(blEnd3 = (blStart3 + blSize3)) 

# Dataframe that reports canonicity for every transcript of the assembly
canonicity <- dplyr::select(assembly,start,end,id) %>%
  mutate(canonicity=ifelse(start>100,"NC", "C")) %>%
  mutate(canonicity=ifelse(end<29000,"NC", canonicity)) %>%
  dplyr::select(-start,-end) 


####################### PROCESSING OF THE DATABASES #######################

################################################################################
############################  PUS7KD vs WT #####################################
################################################################################

# Extraction of the datasets from Nanocompore data
tx_data <- lapply(PUS7_KD_WT_tx, function(x) {
  x <- read_tsv(x, col_types = "ncncccnnncncnc") %>%
    separate(cluster_counts,into = c("SAMPLEID", "Y", "Z"),sep = "_(?=[0-9])",remove = F)%>%
    dplyr::select(-Y,-Z)%>%
    subset(SAMPLEID != "NC")%>%
    mutate(ref_kmer = gsub("T", "U", ref_kmer))%>%
    dplyr::select(-strand,-GMM_cov_type,-GMM_n_clust)
})
  
tx_data <- as.data.frame(bind_rows(tx_data)) 
total_split <- split(tx_data,tx_data$ref_id)
  
toplot <- lapply(X = total_split,FUN = function(x){                             # loop over the transcript models
    x <- x %>%
      mutate(significant=ifelse(GMM_logit_pvalue<=pval_thresh & abs(Logit_LOR)>=LOR_thresh,T,F)) %>% # indicate for every sample using TRUE or FALSE if its LOR and pvalue are significant according to the established thresholds
      rowwise()%>%
      mutate(fleming_presence=get_fleming(genomicPos))%>%   # add a column to indicate if the site is present in the Fleming et al. paper 
      rowwise() %>%
      mutate(ORF = get_ORF(ref_id)) %>%                                           # get the ORF corresponding to the transcript model id
      separate(ref_id, into=c("id","others"), sep="::") %>%
      left_join(canonicity,by="id") %>%                                           # add a column to indicate the canonicity of the transcript model
      rowwise() %>%
      mutate(tx_length = get_tx_length(id)) %>%                                   # assign transcript length
      mutate(start_end_sites = blacklist_start_end(pos,tx_length)) %>%            # assign start and end sites to genomicPos
      left_join(sitelist, by = "ref_kmer") %>%                                    # assign PTM motif to each kmer
      mutate(motif=ifelse(is.na(motif),"Others",motif)) %>%
      rowwise() %>%
      mutate(junction = blacklist_kmer(genomicPos,id)) %>%                        # assign type of junction to each kmer
      rowwise() %>%
      mutate(fragment_ID=assign_fragment_kmer_first_nucl(genomicPos))             # assign RaPID fragment to each kmer
    x$fragment_ID <- factor(x$fragment_ID,levels=c("No_Fragment","Fragment1","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10"))
    return(x)
})

########################### PLOTS #######################################

i<- "PUS7_KD_vs_WT_c37"
pdf(rdp(paste0(i,"_plots_per_transcript.pdf")),height=15,width=20)
lapply(toplot,function(x){
  x %>% 
          {
            ggplot(., aes(x=abs(as.numeric(Logit_LOR)), y=-log10(GMM_logit_pvalue), color=fragment_ID)) +
              geom_point() +
              {if (nrow(subset(x,significant==T))>0) ggrepel::geom_label_repel(data=filter(.,(significant==T)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=5)}+
              scale_color_manual(breaks = c("No_Fragment","Fragment1","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
              ggtitle(unique(x$ORF), subtitle = paste0(unique(x$id)," ",unique(x$canonicity))) +
              theme_bw(22)
          }
})
dev.off()
  
  
pdf(rdp(paste0(i,"_plots_per_transcript_5p.pdf")),height=15,width=20)
lapply(toplot,function(x){
  x %>%
          subset(genomicPos<=100)%>%
          {
            ggplot(., aes(x=abs(as.numeric(Logit_LOR)), y=-log10(GMM_logit_pvalue),color=fragment_ID)) +
              geom_point() +
              {if (nrow(subset(x,significant==T & genomicPos<=100))>0) ggrepel::geom_label_repel(data=filter(., significant==T & genomicPos<=100) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=5)}+
              scale_color_manual(breaks = c("No_Fragment","Fragment1","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
              ggtitle(unique(x$ORF), subtitle = paste0(unique(x$id)," ",unique(x$canonicity))) +
              theme_bw(22)
          }
})
dev.off()

pdf(rdp(paste0(i,"_pvalue_per_position_3p.pdf")),height=15,width=20)
lapply(toplot,function(x){
  x %>%
    subset(genomicPos>=100)%>%
    {
      ggplot(., aes(x=genomicPos, y=-log10(GMM_logit_pvalue),color=fragment_ID)) +
        geom_point() +
        scale_color_manual(breaks = c("No_Fragment","Fragment1","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
        ggtitle(unique(x$ORF), subtitle = paste0(unique(x$id)," ",unique(x$canonicity))) +
        theme_bw(22)
    }
})
dev.off()
############################# PAPER FIGURES ####################################
  
fig_plots <- bind_rows(toplot) %>% subset(canonicity=="C") 
fig_plots <- split(fig_plots,fig_plots$id)
  plots <- lapply(fig_plots,function(x){
    x %>%
      {
        ggplot(., aes(x=abs(as.numeric(Logit_LOR)), y=-log10(GMM_logit_pvalue), color=fragment_ID)) +
          geom_point() +
          {if (nrow(subset(x,significant==T))>0) ggrepel::geom_label_repel(data=filter(.,(significant==T)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=2)}+
          scale_color_manual(breaks = c("No_Fragment","Fragment1","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
          ggtitle(unique(x$ORF), subtitle = paste0(unique(x$id)," ",unique(x$canonicity))) +
          theme_bw()
      }
})
pdf(rdp(paste0(i,"_sharkfin_canonical.pdf")),height=30,width=25)
print(ggarrange(plotlist=plots, ncol = 2,nrow = 3))
dev.off()





############################# WRITE TABLES #####################################


# Table with identity of shared modifications
final_id <- lapply(toplot,function(x){                                          
  x <- x %>%
    subset(significant==T) %>%
    dplyr::select(-SAMPLEID,-start_end_sites,-modif,-others)
  return(x)
})
final_id_all <- as.data.frame(bind_rows(final_id)) %>% 
  mutate(genomicPos=(genomicPos+3))                                             # ALL MODS : significant,no junctions
write.xlsx(final_id_all, rdp(paste0(i,"_modified_sites_",which_nucl,".xls")),sheetName="all_significant_sites",row.names=F,col.names=T)
  
final_id_U <- as.data.frame(bind_rows(final_id)) %>% 
  subset(grepl(which_nucl,ref_kmer)==T)%>%
  subset(canonicity=="C")%>%
  mutate(genomicPos=(genomicPos+3))                                             # ALL CANONICAL Us : significant,no junctions, canonical, U present
write.xlsx(final_id_U, rdp(paste0(i,"_modified_sites_",which_nucl,".xls")),sheetName=paste0("all_canonical_",which_nucl,"s_significant_sites"),row.names=F,col.names=T,append=TRUE)
  
final_id_5p <- as.data.frame(bind_rows(final_id)) %>% 
  subset(genomicPos<=100)%>%
  mutate(genomicPos=(genomicPos+3))                                             # 5p sites : significant,genomicpos <=100
write.xlsx(final_id_5p, rdp(paste0(i,"_modified_sites_",which_nucl,".xls")),sheetName="5p_significant_sites",row.names=F,col.names=T,append=TRUE)
final_id_5p_Us <- as.data.frame(bind_rows(final_id)) %>% 
  subset(genomicPos<=100) %>%
  subset(canonicity=="C") %>%
  subset(grepl(which_nucl,ref_kmer)==T)%>%
  mutate(genomicPos=(genomicPos+3))                                             # 5p Us : significant,U present, canonical, genomicpos <=100
write.xlsx(final_id_5p_Us, rdp(paste0(i,"_modified_sites_",which_nucl,".xls")),sheetName=paste0("5p_canonical_",which_nucl,"s_sign_sites"),row.names=F,col.names=T,append=TRUE)




################################################################################
############################  WT vs IVT #####################################
################################################################################

# Extraction of the datasets from Nanocompore data
tx_data <- lapply(WT_IVT_tx, function(x) {
  x <- read_tsv(x, col_types = "ncncccnnncncnc") %>%
    separate(cluster_counts,into = c("SAMPLEID", "Y", "Z"),sep = "_(?=[0-9])",remove = F)%>%
    dplyr::select(-Y,-Z)%>%
    subset(SAMPLEID != "NC")%>%
    mutate(ref_kmer = gsub("T", "U", ref_kmer))%>%
    dplyr::select(-strand,-GMM_cov_type,-GMM_n_clust)
})

tx_data <- as.data.frame(bind_rows(tx_data)) 
total_split <- split(tx_data,tx_data$ref_id)

toplot <- lapply(X = total_split,FUN = function(x){                             # loop over the transcript models
  x <- x %>%
    mutate(significant=ifelse(GMM_logit_pvalue<=pval_thresh & abs(Logit_LOR)>=LOR_thresh,T,F)) %>% # indicate for every sample using TRUE or FALSE if its LOR and pvalue are significant according to the established thresholds
    rowwise()%>%
    mutate(fleming_presence=get_fleming(genomicPos))%>%   # add a column to indicate if the site is present in the Fleming et al. paper 
    rowwise() %>%
    mutate(ORF = get_ORF(ref_id)) %>%                                           # get the ORF corresponding to the transcript model id
    separate(ref_id, into=c("id","others"), sep="::") %>%
    left_join(canonicity,by="id") %>%                                           # add a column to indicate the canonicity of the transcript model
    rowwise() %>%
    mutate(tx_length = get_tx_length(id)) %>%                                   # assign transcript length
    mutate(start_end_sites = blacklist_start_end(pos,tx_length)) %>%            # assign start and end sites to genomicPos
    left_join(sitelist, by = "ref_kmer") %>%                                    # assign PTM motif to each kmer
    mutate(motif=ifelse(is.na(motif),"Others",motif)) %>%
    rowwise() %>%
    mutate(junction = blacklist_kmer(genomicPos,id)) %>%                        # assign type of junction to each kmer
    rowwise() %>%
    mutate(fragment_ID=assign_fragment_kmer_first_nucl(genomicPos))             # assign RaPID fragment to each kmer
  x$fragment_ID <- factor(x$fragment_ID,levels=c("No_Fragment","Fragment1","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10"))
  return(x)
})

########################### PLOTS #######################################

i<- "WT_IVT"
pdf(rdp(paste0(i,"_plots_per_transcript.pdf")),height=15,width=20)
lapply(toplot,function(x){
  x %>% 
    {
      ggplot(., aes(x=abs(as.numeric(Logit_LOR)), y=-log10(GMM_logit_pvalue), color=fragment_ID)) +
        geom_point() +
        {if (nrow(subset(x,junction=="No junction" & significant==T))>0) ggrepel::geom_label_repel(data=filter(.,(junction=="No junction" & significant==T)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=5)}+
        scale_color_manual(breaks = c("No_Fragment","Fragment1","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
        ggtitle(unique(x$ORF), subtitle = paste0(unique(x$id)," ",unique(x$canonicity))) +
        theme_bw(22)
    }
})
dev.off()


pdf(rdp(paste0(i,"_plots_per_transcript_5p.pdf")),height=15,width=20)
lapply(toplot,function(x){
  x %>%
    subset(genomicPos<=100)%>%
    {
      ggplot(., aes(x=abs(as.numeric(Logit_LOR)), y=-log10(GMM_logit_pvalue),color=fragment_ID)) +
        geom_point() +
        {if (nrow(subset(x,significant==T & genomicPos<=100))>0) ggrepel::geom_label_repel(data=filter(., significant==T & genomicPos<=100) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=5)}+
        scale_color_manual(breaks = c("No_Fragment","Fragment1","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
        ggtitle(unique(x$ORF), subtitle = paste0(unique(x$id)," ",unique(x$canonicity))) +
        theme_bw(22)
    }
})
dev.off()
############################# PAPER FIGURES ####################################

fig_plots <- bind_rows(toplot) %>% subset(canonicity=="C") 
fig_plots <- split(fig_plots,fig_plots$id)
plots <- lapply(fig_plots,function(x){
  x %>%
    {
      ggplot(., aes(x=abs(as.numeric(Logit_LOR)), y=-log10(GMM_logit_pvalue), color=fragment_ID)) +
        geom_point() +
        {if (nrow(subset(x,junction=="No junction" & significant==T))>0) ggrepel::geom_label_repel(data=filter(.,(junction=="No junction" & significant==T)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=2)}+
        scale_color_manual(breaks = c("No_Fragment","Fragment1","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
        ggtitle(unique(x$ORF), subtitle = paste0(unique(x$id)," ",unique(x$canonicity))) +
        theme_bw()
    }
})
pdf(rdp(paste0(i,"_sharkfin_canonical.pdf")),height=30,width=25)
print(ggarrange(plotlist=plots, ncol = 2,nrow = 3))
dev.off()





############################# WRITE TABLES #####################################


# Table with identity of shared modifications
final_id <- lapply(toplot,function(x){                                          
  x <- x %>%
    subset(significant==T) %>%
    dplyr::select(-SAMPLEID,-start_end_sites,-modif,-others)
  return(x)
})
final_id_all <- as.data.frame(bind_rows(final_id)) %>% 
  subset(junction=="No junction")%>%
  mutate(genomicPos=(genomicPos+3))                                             # ALL MODS : significant,no junctions
write.xlsx(final_id_all, rdp(paste0(i,"_modified_sites_",which_nucl,".xls")),sheetName="all_significant_sites",row.names=F,col.names=T)

final_id_U <- as.data.frame(bind_rows(final_id)) %>%
  subset(junction=="No junction")%>%
  subset(grepl(which_nucl,ref_kmer)==T)%>%
  subset(canonicity=="C")%>%
  mutate(genomicPos=(genomicPos+3))                                             # ALL CANONICAL Us : significant,no junctions, canonical, U present
write.xlsx(final_id_U, rdp(paste0(i,"_modified_sites_",which_nucl,".xls")),sheetName=paste0("all_canonical_",which_nucl,"s_significant_sites"),row.names=F,col.names=T,append=TRUE)

final_id_5p <- as.data.frame(bind_rows(final_id)) %>% 
  subset(genomicPos<=100)%>%
  mutate(genomicPos=(genomicPos+3))                                             # 5p sites : significant,genomicpos <=100
write.xlsx(final_id_5p, rdp(paste0(i,"_modified_sites_",which_nucl,".xls")),sheetName="5p_significant_sites",row.names=F,col.names=T,append=TRUE)
final_id_5p_Us <- as.data.frame(bind_rows(final_id)) %>% 
  subset(genomicPos<=100) %>%
  subset(canonicity=="C") %>%
  subset(grepl(which_nucl,ref_kmer)==T)%>%
  mutate(genomicPos=(genomicPos+3))                                             # 5p Us : significant,U present, canonical, genomicpos <=100
write.xlsx(final_id_5p_Us, rdp(paste0(i,"_modified_sites_",which_nucl,".xls")),sheetName=paste0("5p_canonical_",which_nucl,"s_sign_sites"),row.names=F,col.names=T,append=TRUE)

