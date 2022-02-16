library(tidyverse)
library(tidyr)
library(R.utils)
library(reshape2)
library(rlang)
library(stringr)
library(plyr)
library(dplyr)
library(xlsx)



ROOTDIR="/Volumes/scratch/TSSM/cugolini/cov"
RESULTSDIR="/Volumes/scratch/FN/TL/cugolini/cov/scripts/downstream/results_allfiles_LOR05_pval001"
DATADIR="/Volumes/scratch/FN/TL/cugolini/cov/analysis"
DATADIR2="/Volumes/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/nanocompore/comparison"
dir.create(paste0(RESULTSDIR,"/shared/bedtracks"))


########################### PARAMETERS ##########################################

`%!in%` <- Negate(`%in%`)

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

tx <- read_tsv(bdp("analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/orf_annotate/orf_annotate.bed"), col_names=c("chr", "start", "end", "name", "score", "strand", "cdsStart", "cdsEnd", ".", "ex", "exLen", "exSt"), col_types="cnncncnnnncc")
assembly <-read.table(bdp("scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus.bed"), col.names = c("chrom","start","end","id","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"), sep="\t")
file_list <- list.files(path = RESULTSDIR,pattern = "*_modified_sites.xls" , full.names = TRUE,  recursive = F)

# Dataframe that reports canonicity for every transcript of the assembly
canonicity <- dplyr::select(assembly,start,end,id) %>%
  mutate(canonicity=ifelse(start>100,"NC", "C")) %>%
  mutate(canonicity=ifelse(end<29000,"NC", canonicity)) %>%
  mutate(canonicity=ifelse(id=="efad7b96-ac2e-4ce1-9b83-863ffdb18eac|116","NC", canonicity)) %>%  #ORF10
  mutate(canonicity=ifelse(id=="de81ef19-655d-4ced-a9cc-cb8384001058|107","NC", canonicity)) %>%  #ORF9D
  dplyr::select(-start,-end) 

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


assembly <- left_join(assembly,canonicity,by="id")
names <- names %>% select(id,name) %>% separate(id,into=c("id"),sep = "::",remove=T)
assembly <- left_join(assembly,names,by="id")

assembly_multiple_isoform <- assembly %>% 
  subset(canonicity=="C") %>% 
  group_by(name) %>% 
  filter(n()>1)

####################### PROCESSING OF THE DATABASES ############################

sites_all <- lapply(file_list,function(x){
  linename <- str_remove(basename(x),"_modified_sites.xls")
  x <- read.xlsx(x,sheetName ="all_significant_sites",as.data.frame = T)
  x <- x %>%
    select(pos,chr,genomicPos,ref_id,others,ref_kmer,ORF,canonicity,`Modification.type`,IVT,fragment_ID,burrows_presence)%>%
    dplyr::rename(Motif=Modification.type)%>%
    mutate(cell_type=linename)
})
cols_selection <- c("pos","chr","genomicPos","ref_id","others","ref_kmer","ORF","canonicity","Motif","IVT","fragment_ID","burrows_presence")

sites_all_1_isoform <- lapply(sites_all, function(x){x%>%select(-cell_type)})
sites_all_1_isoform <-plyr::join_all(sites_all_1_isoform, by=cols_selection,type="inner") %>%
  subset(ref_id %!in% assembly_multiple_isoform$id)

sites_all_multiple_isoform <- lapply(sites_all, function(x){
  x<- x %>%
    subset(ref_id %in% assembly_multiple_isoform$id) 
})
sites_all_multiple_isoform <- bind_rows(sites_all_multiple_isoform)
sites_all_multiple_isoform <- split(sites_all_multiple_isoform,sites_all_multiple_isoform$genomicPos)
sites_all_multiple_isoform <- lapply(sites_all_multiple_isoform,function(x){
   x <- split(x,x$ORF)
   x <- lapply(x,function(y){
     if(length(unique(y$cell_type))==3){y<-y%>%slice(1:1)%>%mutate(del=F)}
     else{y<-y%>%mutate(del=T)}
  })
   x <- bind_rows(x)
})
sites_all_multiple_isoform <- bind_rows(sites_all_multiple_isoform) %>%
  subset(del==F)%>%
  select(chr,genomicPos,ref_kmer,ORF,canonicity,Motif,IVT,fragment_ID,burrows_presence)%>%
  mutate(pos="NA")%>%
  mutate(ref_id="NA")%>%
  mutate(others="NA")
 
sites_all<- rbind(sites_all_1_isoform,sites_all_multiple_isoform)               # all significant sites                                
canonical_signif_U <- sites_all %>%                                             # canonical U sites (no junction)
  subset(IVT=="No junction") %>% 
  subset(substring(ref_kmer,3,3)=="U")%>%
  subset(canonicity=="C")
signif_5p <- sites_all %>%                                                      # significant sites before 100 
  subset((genomicPos-3)<=100)
signif_5p_canonical_U <- sites_all %>%                                          # significant canonical U sites before 100 
  subset((genomicPos-3)<=100) %>%
  subset(canonicity=="C") %>%
  subset(substring(ref_kmer,3,3)=="U")

### peak called datasets

file_list <- list.files(path = RESULTSDIR,pattern = "*_modified_sites.xls" , full.names = TRUE,  recursive = F)
peakcalled <- lapply(file_list,function(x){
  x <- read.xlsx(x,sheetName ="peakcalled_canonical_Us",as.data.frame = T)
})
cols_selection <- c("pos","chr","genomicPos","ref_id","strand","ORF","canonicity","IVT","fragment_ID","burrows_presence")

peakcalled_1_isoform <- peakcalled %>%
  bind_rows() %>%
  subset(ref_id %!in% assembly_multiple_isoform$id) %>%
  distinct()
appo <- lapply(peakcalled, function(x){ x %>%
    subset(ref_id %!in% assembly_multiple_isoform$id)})
present_in_caco2 <- do.call(paste0, peakcalled_1_isoform) %in% do.call(paste0, appo[[1]])
present_in_calu3 <- do.call(paste0, peakcalled_1_isoform) %in% do.call(paste0, appo[[2]])
present_in_vero <- do.call(paste0, peakcalled_1_isoform) %in% do.call(paste0, appo[[3]])

peakcalled_1_isoform <- cbind(peakcalled_1_isoform,present_in_caco2,present_in_calu3,present_in_vero)


peakcalled <- lapply(file_list,function(x){
  x <- read.xlsx(x,sheetName ="peakcalled_canonical_Us",as.data.frame = T) %>%
    subset(ref_id %in% assembly_multiple_isoform$id)
})
peakcalled <- lapply(peakcalled,function(x){
  x <- split(x, x$ORF)
  x <- lapply(x, function(y){
    appo <- split(y,y$ref_id) %>% lapply(function(z){z<-z%>%select(-pos,-ref_id)})
    y <- y %>% select(-pos,-ref_id)
    iso1 <- do.call(paste0, y) %in% do.call(paste0, appo[[1]])
    iso2 <- do.call(paste0, y) %in% do.call(paste0, appo[[2]])
    iso3 <- do.call(paste0, y) %in% do.call(paste0, appo[[3]])
    y <- cbind(y,iso1,iso2,iso3) %>% 
      subset(iso1!=FALSE | iso2!=FALSE | iso3!=FALSE) %>% 
      select(-iso1,-iso2,-iso3)
    return(y)
  })
  x <- bind_rows(x)
})

appo <- peakcalled
peakcalled_multiple_isoform <- peakcalled %>%
  bind_rows() %>%
  distinct()
present_in_caco2 <- do.call(paste0, peakcalled_multiple_isoform) %in% do.call(paste0, appo[[1]])
present_in_calu3 <- do.call(paste0, peakcalled_multiple_isoform) %in% do.call(paste0, appo[[2]])
present_in_vero <- do.call(paste0, peakcalled_multiple_isoform) %in% do.call(paste0, appo[[3]])

peakcalled_multiple_isoform <- cbind(peakcalled_multiple_isoform,present_in_caco2,present_in_calu3,present_in_vero) %>%
  mutate(pos="NA")%>%
  mutate(ref_id="NA")
final <- rbind(peakcalled_1_isoform,peakcalled_multiple_isoform) %>%
  as.data.frame()%>%
  rowwise()%>%
  mutate(presence=sum(c(present_in_caco2,present_in_calu3,present_in_vero), na.rm = TRUE))%>%
  as.data.frame()

# Burrows

burrows_redundant_signif <- sites_all %>%                                       # all Burrows significant sites
  subset(burrows_presence==T) 
burrows_nonredundant_signif <- sites_all %>%                                  # all Burrows significant non redundant sites
  subset(burrows_presence==T) %>%
  select(genomicPos,ref_kmer)%>%
  dplyr::group_by(genomicPos,ref_kmer) %>% 
  dplyr::filter(row_number() == 1)


burrows_all <- lapply(file_list,function(x){
  linename <- str_remove(basename(x),"_modified_sites.xls")
  x <- read.xlsx(x,sheetName ="Burrows_redundant_ALL_sites",as.data.frame = T)
  x <- x %>%
    select(pos,chr,genomicPos,ref_id,others,ref_kmer,ORF,canonicity,`Modification.type`,IVT,fragment_ID,burrows_presence)%>%
    dplyr::rename(Motif=Modification.type)%>%
    mutate(cell_type=linename)
})
cols_selection <- c("pos","chr","genomicPos","ref_id","others","ref_kmer","ORF","canonicity","Motif","IVT","fragment_ID","burrows_presence")

burrows_all_1_isoform <- lapply(burrows_all, function(x){x%>%select(-cell_type)})
burrows_all_1_isoform <-plyr::join_all(burrows_all_1_isoform, by=cols_selection,type="inner") %>%
  subset(ref_id %!in% assembly_multiple_isoform$id)

burrows_all_multiple_isoform <- lapply(burrows_all, function(x){
  x<- x %>%
    subset(ref_id %in% assembly_multiple_isoform$id) 
})
burrows_all_multiple_isoform <- bind_rows(burrows_all_multiple_isoform)
burrows_all_multiple_isoform <- split(burrows_all_multiple_isoform,burrows_all_multiple_isoform$genomicPos)
burrows_all_multiple_isoform <- lapply(burrows_all_multiple_isoform,function(x){
  x <-split(x,x$ORF)
  x <- lapply(x,function(y){
    if(length(unique(y$cell_type))==3){y<-y%>%slice(1:1)%>%mutate(del=F)}
    else{y<-y%>%mutate(del=T)}
  })
  x <- bind_rows(x)
})
burrows_all_multiple_isoform <-bind_rows(burrows_all_multiple_isoform) %>%
  subset(del==F)%>%
  select(chr,genomicPos,ref_kmer,ORF,canonicity,Motif,IVT,fragment_ID,burrows_presence)%>%
  mutate(pos="NA")%>%
  mutate(ref_id="NA")%>%
  mutate(others="NA")

burrows_all<- rbind(burrows_all_1_isoform,burrows_all_multiple_isoform)         # all Burrows sites
burrows_nonredundant_all <- burrows_all %>%                                     # all Burrows non redundant sites
  dplyr::select(genomicPos,ref_kmer)%>%
  dplyr::group_by(genomicPos,ref_kmer) %>% 
  dplyr::filter(row_number() == 1)
burrows_nonredundant_all <- as.data.frame(burrows_nonredundant_all)





############################# WRITE TABLES #####################################

# Table with identity of shared modifications

write.xlsx(sites_all, rdp("shared_sites.xls"),sheetName="all_significant_sites",row.names=F,col.names=T)
write.xlsx(canonical_signif_U, rdp("shared_sites.xls"),sheetName="all_canonical_Us_significant_sites",row.names=F,col.names=T,append=TRUE)
write.xlsx(signif_5p, rdp("shared_sites.xls"),sheetName="5p_significant_sites",row.names=F,col.names=T,append=TRUE)
write.xlsx(signif_5p_canonical_U, rdp("shared_sites.xls"),sheetName="5p_canonical_Us_significant_sites",row.names=F,col.names=T,append=TRUE)
write.xlsx(final, rdp("shared_sites.xls"),sheetName="peakcalled_Us",row.names=F,col.names=T,append=TRUE)

# Table with Burrows sites
write.xlsx(burrows_all, rdp("shared_sites.xls"),sheetName="Burrows_redundant_ALL_sites",row.names=F,col.names=T,append=TRUE)
write.xlsx(burrows_redundant_signif, rdp("shared_sites.xls"),sheetName="Burrows_redundant_sign_sites",row.names=F,col.names=T,append=TRUE)
write.xlsx(burrows_nonredundant_signif, rdp("shared_sites.xls"),sheetName="Burrows_non_redundant_CandNC_sign_sites",row.names=F,col.names=T,append=TRUE)
write.xlsx(burrows_nonredundant_all, rdp("shared_sites.xls"),sheetName="Burrows_non_redundant_CandNC_all_sites",row.names=F,col.names=T,append=TRUE)



###  Working on high-confidence pseudoU sites

f <- list.files(rdp("shared/"), include.dirs = F, full.names = T, recursive = T)# remove the files
file.remove(f)
canonical_signif_U <- split(canonical_signif_U,canonical_signif_U$genomicPos)
lapply(canonical_signif_U,function(x){
  genPos <- unique(x$genomicPos)
  ref_kmer <- unique(x$ref_kmer)
  x <- x %>%
    mutate(ref_tx=paste(ref_id,"::",others,sep = ""))%>%
    select(pos,ref_tx)
  write.table(x,sep = "\t",quote = F,row.names = F,col.names = F,file = rdp(paste("shared/",ref_kmer,"_",genPos,".txt",sep="")))
})


### Build genomic track


final2 <- final %>% 
  subset(presence>=2)
final3 <- final %>% 
  subset(presence==3)
final2 <- split(final2,final2$ORF)
final3 <- split(final3,final3$ORF)
bed_list <- lapply(final2, function(x){
  ORF_id <- unique(x$ORF)
  x <- x %>%
    mutate(start=genomicPos) %>%
    mutate(end=genomicPos) %>%
    mutate(name=paste0("U_genomicpos=",genomicPos,"_ORF=",ORF_id))%>%
    mutate(score=0)%>%
    mutate(strand="+") %>%
    select(chr,start,end,name,score,strand)
  write.table(x,sep = "\t",quote = F,row.names = F,col.names = F,file = rdp(paste("shared/bedtracks/2_shared_peakcalled_Us_",ORF_id,".bed",sep="")))
  
})

bed_list <- lapply(final3, function(x){
  ORF_id <- unique(x$ORF)
  x <- x %>%
    mutate(start=genomicPos) %>%
    mutate(end=genomicPos) %>%
    mutate(name=paste0("U_genomicpos=",genomicPos,"_ORF=",ORF_id))%>%
    mutate(score=0)%>%
    mutate(strand="+") %>%
    select(chr,start,end,name,score,strand)
  write.table(x,sep = "\t",quote = F,row.names = F,col.names = F,file = rdp(paste("shared/bedtracks/3_shared_peakcalled_Us_",ORF_id,".bed",sep="")))
  
})


# Canonical assembly track
assembly_track <- assembly %>%
  subset(canonicity=="C") %>%
  subset(id %!in% assembly_multiple_isoform$id)%>%
  select(-canonicity,-id)%>%
  relocate(name, .after = end)
assembly_multiple_isoform <- head(assembly_multiple_isoform,2) %>%
  select(-canonicity,-id)%>%
  relocate(name, .after = end)
assembly_track <- rbind(assembly_track,assembly_multiple_isoform)
write.table(assembly_track,sep = "\t",quote = F,row.names = F,col.names = F,file = rdp("shared/bedtracks/assembly.bed"))

