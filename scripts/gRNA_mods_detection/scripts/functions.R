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

# Function that returns full path from filesdir
fdp <- function(relpath){
  return(paste0(FILES,"/",relpath))
}

# Function that returns full path from currentdir
cdp <- function(relpath){
  return(paste0(CURRDIR,"/",relpath))
}

# Functions for plot
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

vlab <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
cbbPalette <- c("#0072B2","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7")


# Function to know if Fleming et al. sites are present
get_fleming = function(kmer_pos) {
  kmer <- seq(kmer_pos,kmer_pos+4,1)
  if(is_empty(intersect(kmer,fleming_sites_bed)) != T) {return(T)
  }else{return(F)}
}


# Function to get overlap of a sequence interval with fragments (even if partial)
get_fragment_partial_overlap = function(kmer_pos,fragments_coordfile){
  kmer <- seq(kmer_pos,kmer_pos+4,1)
  fragments_coordfile <- fragments_coordfile %>%
    rowwise() %>%
    mutate(overlap=length(intersect(kmer,seq(start,end,1))))
  if(max(fragments_coordfile$overlap) == 0){
    return("No_Fragment")
  }else{
    return(fragments_coordfile[which.max(fragments_coordfile$overlap),]$ID)
  }
}


# obtain a vector of nucleotides overlapping extremities of a bed file according to specific intervals
extr_overlap = function(coordfile,left_interval,right_interval){
  junctions <- coordfile %>% 
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
kmer_overlaps_IVT = function(kmer_pos,IVT_coord_file,left_int,right_int){
  kmer <- seq(kmer_pos,kmer_pos+4,1)
  nucl_overlap_IVT_junction <- extr_overlap(IVT_coord_file,left_int,right_int)
  if(is_empty(intersect(kmer,nucl_overlap_IVT_junction)) != T) {return("IVT_junction")
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
  int_vector <- apply(sgRNAs_annotation, 1, function(x) length(intersect(seq(x['bed_start'],(as.numeric(x['bed_end'])-1),1),kmer)))
  if(max(int_vector)==0){return(F)}else{return(T)}
}

# Negate operator "in"
`%!in%` <- Negate(`%in%`)
