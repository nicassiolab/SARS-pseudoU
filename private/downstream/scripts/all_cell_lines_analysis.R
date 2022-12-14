library(tidyverse)
library(ggpubr)
library(GGally)
library(tidyr)
library(R.utils)
library(reshape2)
library(rlang)
library(stringr)
library(plyr)
library(dplyr)



ROOTDIR="/Volumes/scratch/TSSM/cugolini/cov"
RESULTSDIR="/Volumes/scratch/FN/TL/cugolini/cov/scripts/downstream/results"
DATADIR="/Volumes/scratch/FN/TL/cugolini/cov/analysis"
DATADIR2="/Volumes/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/nanocompore/comparison"
dir.create(paste0(RESULTSDIR,"/shared"))

########################### PARAMETERS ##########################################



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


########################### GENERAL DATA #######################################

sites_all <- list.files(path = RESULTSDIR,pattern = "*_sites_identity.txt" , full.names = TRUE,  recursive = F)
sites_5p <- list.files(path = RESULTSDIR,pattern = "*_sites_identity_5p.txt" , full.names = TRUE,  recursive = F)
sites_burrows <- list.files(path = RESULTSDIR,pattern = "*_sites_identity_burrows.txt" , full.names = TRUE,  recursive = F)

####################### PROCESSING OF THE DATABASES ############################

newnames <- c("final_pvalue","final_LOR")

sites_all <- lapply(sites_all,function(x){
  linename <- str_remove(basename(x),"_sites_identity.txt")
  x <- read_tsv(x)
  x <- x %>%
    select(pos,chr,genomicPos,ref_id,others,ref_kmer,ORF,canonicity,`Modification type`,IVT,fragment_ID,final_pvalue,final_LOR,burrows_presence)%>%
    mutate(on_the_U=ifelse(substr(ref_kmer,3,3)=="U",T,F))%>%
    subset(on_the_U==T)%>%
    rename_with(~ paste0(.x,"_",linename),.cols = newnames)
})

cols_selection <- c("pos","chr","genomicPos","ref_id","others","ref_kmer","ORF","canonicity","Modification type","IVT","fragment_ID","burrows_presence","on_the_U")
sites_all <- join_all(sites_all, by=cols_selection,type="inner")
write.table(sites_all, rdp(paste0("/shared/","sites_identity.txt")),sep="\t",quote=F,row.names=F,col.names=T)


sites_5p <- lapply(sites_5p,function(x){
  linename <- str_remove(basename(x),"_sites_identity_5p.txt")
  x <- read_tsv(x)
  x <- x %>%
    select(pos,chr,genomicPos,ref_id,others,ref_kmer,ORF,canonicity,`Modification type`,IVT,fragment_ID,final_pvalue,final_LOR,burrows_presence)%>%
    mutate(on_the_U=ifelse(substr(ref_kmer,3,3)=="U",T,F))%>%
    subset(on_the_U==T)%>%
    rename_with(~ paste0(.x,"_",linename),.cols = newnames)
})

cols_selection <- c("pos","chr","genomicPos","ref_id","others","ref_kmer","ORF","canonicity","Modification type","IVT","fragment_ID","burrows_presence","on_the_U")
sites_5p <- join_all(sites_5p, by=cols_selection,type="inner")
write.table(sites_5p, rdp(paste0("/shared/","sites_identity_5p.txt")),sep="\t",quote=F,row.names=F,col.names=T)


newnames <- c("final_pvalue","final_LOR","shared")
sites_burrows <- lapply(sites_burrows,function(x){
  linename <- str_remove(basename(x),"_burrows_sites_identity.txt")
  x <- read_tsv(x)
  x <- x %>%
    select(pos,chr,genomicPos,ref_id,others,ref_kmer,ORF,canonicity,`Modification type`,IVT,fragment_ID,final_pvalue,final_LOR,burrows_presence,shared)%>%
    mutate(on_the_U=ifelse(substr(ref_kmer,3,3)=="U",T,F))%>%
    subset(on_the_U==T)%>%
    rename_with(~ paste0(.x,"_",linename),.cols = newnames)
})

cols_selection <- c("pos","chr","genomicPos","ref_id","others","ref_kmer","ORF","canonicity","Modification type","IVT","fragment_ID","burrows_presence","on_the_U")
sites_burrows <- join_all(sites_burrows, by=cols_selection,type="inner")
write.table(sites_burrows, rdp(paste0("/shared/","burrows_sites_identity.txt")),sep="\t",quote=F,row.names=F,col.names=T)


