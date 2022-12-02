library(purrr)
library(DescTools)

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
ddp2 <- function(relpath){
  return(paste0(DATADIR2,"/",relpath))
}

vlab <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

# Function to get the PTM motif in a genomic sequence
get_motif = function(which_nucl,left_int,right_int,motifs_df,genome_fa){
  if(which_nucl=="U"){
    motifs <- vector()
    for(i in 1:nrow(motifs_df)){
      if(motifs_df[i,"motif"]=="UA"){
        seq_to_search <- gsub("T","U",toupper(substr(genome_fa,left_int,right_int)))  # from DNA to RNA and exclude extremities (for the modified U)
        element <- unlist(str_extract_all(seq_to_search,as.character(motifs_df[i,"ref_kmer"]))) # search ref_kmer over the fasta
        if(is_empty(element)==F){motifs <-append(motifs,element)}
      }else if(motifs_df[i,"motif"]=="UNUAR"){
        seq_to_search <- gsub("T","U",toupper(substr(genome_fa,(left_int-2),(right_int+2))))
        element <- unlist(str_extract_all(seq_to_search,as.character(motifs_df[i,"ref_kmer"])))
        if(is_empty(element)==F){motifs <- append(motifs,element)}
      }else if(motifs_df[i,"motif"]=="GUUNNA"){
        seq_to_search <- gsub("T","U",toupper(substr(genome_fa,(left_int-2),(right_int+3))))
        element <- unlist(str_extract_all(seq_to_search,as.character(motifs_df[i,"ref_kmer"])))
        if(is_empty(element)==F){motifs <-append(motifs,element)}
      }
    }
    return(toString(motifs))
  }else if(which_nucl=="A"){
      motifs <- vector()
      for(i in 1:nrow(motifs_df)){
        seq_to_search <- gsub("T","U",toupper(substr(genome_fa,left_int-2,right_int+2)))
        element <- unlist(str_extract_all(seq_to_search,as.character(motifs_df[i,"ref_kmer"])))
        if(is_empty(element)==F){motifs <-append(motifs,element)}
      }
    return(toString(motifs))
  }
}


# Function to know if Fleming et al. sites are present
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
get_fragment_partial_overlap = function(left_int,right_int,fragments_bedfile){
  fragments_bedfile <- fragments_bedfile %>%
    mutate(start=as.numeric(start),end=as.numeric(end)) %>%
    rowwise %>%
    mutate(overlap=length(intersect(seq(left_int,right_int,1),seq(start,end,1))))
  if(max(fragments_bedfile$overlap) == 0){
    return("No_Fragment")
  }else{
    return(fragments_bedfile[which.max(fragments_bedfile$overlap),]$ID)
  }
}


# obtain a vector of nucleotides overlapping extremities of a bed file according to specific intervals
extr_overlap = function(bedfile,left_interval, right_interval){
  junctions <- bedfile %>% 
    dplyr::select(start,end) %>% 
    unlist(use.names=FALSE) %>% 
    sort()
  blacklist <- vector()
  for (i in junctions){
    temp <- seq(from = i-left_interval, to = i+right_interval, by = 1)
    blacklist <- c(blacklist,temp)
  }
  blacklist <- blacklist[which(blacklist >= 0)]
  return(blacklist)
}

# check if any of the nucleotides composing the k-mer overlaps an IVT junction
kmer_overlaps_IVT = function(kmer_pos,IVT_bed_file,left_int,right_int){
  kmer <- seq(kmer_pos,kmer_pos+4,1)
  nucl_overlap_IVT_junction <- extr_overlap(IVT_bed_file,left_int,right_int)
  if(is_empty(intersect(kmer,nucl_overlap_IVT_junction)) != T) {return("IVT_junction")
    }else{return("No_junction")}
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

# check if any of the nucleotides composing the k-mer overlaps an IVT or ORF junction
kmer_overlaps_ORFjunc = function(kmer_pos,tx_id){
  kmer <- seq(kmer_pos,kmer_pos+4,1)
  if(is_empty(intersect(kmer,junc_blacklist(tx_id))) != T) {return("ORF_junction")
    }else{return("No_junction")}
}


# Function to annotate the SNPs from different strains
kmer_overlaps_SNPs = function(kmer_pos,SNPs_annotation){
  kmer <- seq(kmer_pos,kmer_pos+4,1)
  if(is_empty(intersect(kmer,as.vector(SNPs_annotation$genomicPos))) != T) {return("SNP")
  }else{return(NULL)}
}


# Function to intersect kmers and sgRNA sites intervals
kmer_overlaps_sgRNA_site = function(kmer_pos,sgRNAs_annotation){
  kmer <- seq(kmer_pos,kmer_pos+4,1)
  int_vector <-apply(sgrna_sites, 1, function(x) length(intersect(seq(x['start'],x['end'],1),kmer)))
  if(max(int_vector)==0){return(F)}else{return(T)}
}

