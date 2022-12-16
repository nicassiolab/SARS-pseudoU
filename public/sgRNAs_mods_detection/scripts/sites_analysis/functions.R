library(tidyverse)
library(tidyr)
library(R.utils)
library(reshape2)
library(rlang)
library(stringr)
library(dplyr)
library(xlsx)
library(seqinr)
library(ggpubr)
library(DescTools)

########################### FUNCTIONS ##########################################
# Function to negate operator in
`%!in%` <- Negate(`%in%`)

# Function that returns full path from basedir
bdp <- function(relpath){
  return(paste0(ROOTDIR,"/",relpath))
}

# Function that returns full path from basedir
fdp <- function(relpath){
  return(paste0(FILES,"/",relpath))
}


# Function that returns full path from resultsdir
rdp <- function(relpath){
  return(paste0(RESULTSDIR,"/",relpath))
}

# Function that returns full path from currentdir
cdp <- function(relpath){
  return(paste0(CURRDIR,"/",relpath))
}



# Function to get the PTM motif in a sequence
if(which_nucl=="U"){
  get_motif = function(left_int,right_int){
    motifs <- vector()
    for(i in 1:nrow(sitelist)){
      if(sitelist[i,"motif"]=="UA"){
        seq_to_search <- gsub("T","U",toupper(substr(viral_genome,left_int,right_int)))
        element <- unlist(str_extract_all(seq_to_search,as.character(sitelist[i,"ref_kmer"])))
        if(is_empty(element)==F){motifs <-append(motifs,element)}
      }else if(sitelist[i,"motif"]=="UNUAR"){
        seq_to_search <- gsub("T","U",toupper(substr(viral_genome,(left_int-2),(right_int+2))))
        element <- unlist(str_extract_all(seq_to_search,as.character(sitelist[i,"ref_kmer"])))
        if(is_empty(element)==F){motifs <- append(motifs,element)}
      }else if(sitelist[i,"motif"]=="GUUNNA"){
        seq_to_search <- gsub("T","U",toupper(substr(viral_genome,(left_int-2),(right_int+3))))
        element <- unlist(str_extract_all(seq_to_search,as.character(sitelist[i,"ref_kmer"])))
        if(is_empty(element)==F){motifs <-append(motifs,element)}
      }
    }
    return(toString(motifs))
  }
}else if(which_nucl=="A"){
  get_motif = function(left_int,right_int){
    motifs <- vector()
    for(i in 1:nrow(sitelist)){
      seq_to_search <- gsub("T","U",toupper(substr(viral_genome,left_int-2,right_int+2)))
      element <- unlist(str_extract_all(seq_to_search,as.character(sitelist[i,"ref_kmer"])))
      if(is_empty(element)==F){motifs <-append(motifs,element)}
    }
    return(toString(motifs))
  }
}

# Function to know if fleming sites are present
get_fleming = function(left_x,right_x) {
  presence <- sapply(fleming_paper_notation, function(x){
    ifelse((x>=left_x & x<=right_x),T,F)
  })
  return(any(presence))
}

# Function to get the canonical tx_id based on the ORF
get_tx = function(ORF) {
  names_C <- subset(assembly, canonicity=="C") %>% select(id,name)
  temp <- names_C[with(names_C, name %in% ORF),]$id
  return(temp)
}

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
get_fleming_single_site = function(site) {
  blacklist_fleming <- vector()
  for (i in fleming_sites_bed){
    temp <- seq(from = i, to = i+4, by = 1)
    blacklist_fleming <- c(blacklist_fleming,temp)
  }
  return(site %in% blacklist_fleming)
}

# Function to get a transcript length
get_tx_length = function(transcript) {
  tx_length = tx_lengths[with(tx_lengths, tx_lengths$id %in% transcript),]$sumrow
  return(tx_length)
}

# Function to get the ORF encoded in a transcript
get_ORF = function(transcript) {
  temp <- names[with(names, id %in% transcript),]$ORF
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






# Function to get the ORF encoded nucleotide per nucleotide
get_ORF_per_nucl = function(ORF_name,left_int,right_int) {
  assembly_C <- subset(assembly,canonicity=="C")
  temp <- assembly_C[with(assembly_C, ORF %in% ORF_name),]$id
  blacklist_tx <- vector()
  for (i in temp){
    blacklist_tx <- c(blacklist_tx,junc_blacklist(i))
  }
  blacklist_tx <- unique(sort(blacklist_tx))
  ORF_junctions <- intersect(blacklist_tx,seq(left_int,right_int,by=1))
  IVT_junctions <- intersect(blacklist_IVT,seq(left_int,right_int,by=1))
  tot_junctions <- unique(sort(c(ORF_junctions,IVT_junctions)))
  return(toString(tot_junctions))
}

# Function to get shared intervals of nucleotides
get_shared = function(data_list) {
  ORF_id<-unique(data_list$ORF)
  full_seq <- seq(min(data_list$genomicPos),max(data_list$genomicPos),1) %>%
    as.data.frame()
  colnames(full_seq)<- c("genomicPos")
  data_list <- as.data.frame(left_join(full_seq,data_list,by="genomicPos"))
  right <- data.frame()
  if(is.na(data_list[2,"ORF"])==T){right <- as.data.frame(data_list[1,"genomicPos"])}
  left <- as.data.frame(data_list[1,"genomicPos"])
  for (row in 2:nrow(data_list)){
    if(is.na(data_list[(row+1),"ORF"])==T & is.na(data_list[(row-1),"ORF"])==F & is.na(data_list[(row),"ORF"])==F){
      right <- rbind(right,data_list[row,"genomicPos"])}
    else if(is.na(data_list[(row+1),"ORF"])==T & is.na(data_list[(row-1),"ORF"])==T & is.na(data_list[(row),"ORF"])==F){
      right <- rbind(right,data_list[row,"genomicPos"])}
  }
  colnames(right) <- c("right_interval")
  colnames(left) <- c("left_interval")
  for (row in 2:nrow(data_list)){
    if(is.na(data_list[(row+1),"ORF"])==F & is.na(data_list[(row-1),"ORF"])==T & is.na(data_list[(row),"ORF"])==F){
      left <- rbind(left,data_list[row,"genomicPos"])}
    else if(is.na(data_list[(row+1),"ORF"])==T & is.na(data_list[(row-1),"ORF"])==T & is.na(data_list[(row),"ORF"])==F){
      left <- rbind(left,data_list[row,"genomicPos"])}
  }
  colnames(left) <- c("left_interval")
  interval <- cbind(left,right)%>%
    rowwise()%>%
    dplyr::mutate(sequence=gsub("T","U",toupper(substr(viral_genome,left_interval, right_interval)))) %>%
    rowwise()%>%
    dplyr::mutate(mod=get_motif(left_interval,right_interval)) %>%
    dplyr::mutate(ORF=ORF_id)
  return(interval)
}

# Function to get presence in multiple cell lines
get_presence = function(peaked_data){
  appo <- peaked_data
  peaked_data <- bind_rows(peaked_data) %>%
    distinct()
  present_in_caco2 <- do.call(paste0, peaked_data) %in% do.call(paste0, appo[[1]])
  present_in_calu3 <- do.call(paste0, peaked_data) %in% do.call(paste0, appo[[2]])
  present_in_vero <- do.call(paste0, peaked_data) %in% do.call(paste0, appo[[3]])
  peaked_data <- cbind(peaked_data,present_in_caco2,present_in_calu3,present_in_vero)%>%
    distinct() %>%
    rowwise() %>%
    dplyr::mutate(presence=sum(c(present_in_caco2,present_in_calu3,present_in_vero), na.rm = TRUE))
}

# Function to get overlap of a sequence interval with fragments (even if partial)
get_fragment_partial_overlap = function(left_int,right_int){
  fragments_vec <- apply(fragments, 1, function(x){
    if(exists("fragments_list")==F){fragments_list <- vector()}
    if((c(left_int,right_int) %overlaps% c(as.numeric(x[2]),as.numeric(x[3])))==T){
      fragments_list <- append(fragments_list,as.character(x[1]))
    }
    return(fragments_list)
  })
  return(toString(unlist(fragments_vec)))
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
