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



ROOTDIR="/Volumes/scratch/TSSM/cugolini/cov"
RESULTSDIR="/Volumes/scratch/FN/TL/cugolini/cov/scripts/downstream/results_allfiles_LOR05_pval001"
DATADIR="/Volumes/scratch/FN/TL/cugolini/cov/analysis"
DATADIR2="/Volumes/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/nanocompore/comparison"
dir.create(paste0(RESULTSDIR,"/shared/bedtracks"),showWarnings = F)
dir.create(paste0(RESULTSDIR,"/shared/corresp_tables"),showWarnings = F)
dir.create(paste0(RESULTSDIR,"/shared/tracks_for_overlap"),showWarnings = F)


########################### PARAMETERS #########################################

`%!in%` <- Negate(`%in%`)
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



########################### GENERAL DATA #######################################

IVT <- read_tsv(bdp("scripts_new/backupped_data/IVT_junctions.bed"), col_names = c("chrom","start","end","id","score","strand"),col_types="cnncnc") %>% 
  dplyr::mutate(start=(start-1),end=(end-1))
assembly <- read.table(bdp("scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus.bed"), col.names = c("chrom","start","end","id","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"), sep="\t")
tx <- read_tsv(bdp("analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/orf_annotate/orf_annotate.bed"), col_names=c("chr", "start", "end", "ref_id", "score", "strand", "cdsStart", "cdsEnd", ".", "ex", "exLen", "exSt"), col_types="cnncncnnnncc")
tx_lengths <- read.table(bdp("analysis/recappable_assembly/two_datasets/assemblies/pinfish/aln_consensus.bed"), sep = '\t',header = FALSE) %>%
  separate(V11, into=c("ex1","ex2","ex3") ,sep=",", remove=F) 
fragments <- read_tsv(bdp("scripts_new/backupped_data/RAPID/fragments_genomic_coord_UCSC.txt"),col_types = "cnn")
fragments_bed <- fragments %>% 
  dplyr::mutate(start=(start-1),end=(end-1))
viral_genome <- read.fasta(file="/Volumes/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa") %>%
  getSequence(as.string=T) %>% unlist()
file_list <- list.files(path = RESULTSDIR,pattern = paste0("*_modified_sites_",which_nucl,".xls") , full.names = TRUE,  recursive = F)


##  input files that are different depending in which nucleotide the modification is
if(which_nucl=="U"){
  sitelist <- read_tsv(bdp("scripts_new/backupped_data/sites_pseudoU_motifs_blacklist.txt"),col_types="ccc") %>%
    dplyr::rename(motif=IUPAC) %>%
    subset(motif=="UNUAR")                                                      # use only UNUAR motifs
}else if(which_nucl=="A"){
  sitelist <- read_tsv(bdp("scripts_new/backupped_data/sites_blacklist.txt"),col_types="ccc") %>%
    dplyr::rename(motif=IUPAC) %>% subset(modif=="m6A")
}

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


assembly <- left_join(assembly,canonicity,by="id")
names <- names %>% select(ref_id,ORF) %>% separate(ref_id,into=c("id","others"),sep = "::",remove=T)
assembly <- left_join(assembly,names,by="id")

assembly_multiple_isoform <- assembly %>% 
  subset(canonicity=="C") %>% 
  group_by(ORF) %>% 
  filter(n()>1)


####################### PROCESSING OF THE DATABASES ############################

### peak called datasets

peakcalled_1_iso <- lapply(file_list,function(x){
  linename <- str_remove(basename(x),paste0("_modified_sites_",which_nucl,".xls"))
  x <- read.xlsx(x,sheetName = paste0("peakcalled_canonical_",which_nucl,"s"),as.data.frame = T)%>%
    mutate(cell_type=linename) %>%
    subset(id %!in% assembly_multiple_isoform$id) %>%
    select(-fleming_presence)
  x <- split(x, x$id)
  x <- lapply(x, function(y){
    tx <- unique(y$id)
    y <-y %>% 
      select("genomicPos","id")%>%
      mutate(type="central")
    tot_new <- data.frame()
    for (row in 1:nrow(y)){
      genpos <- y[row, "genomicPos"]
      new_sites <- as.data.frame(seq((genpos-4),(genpos+4),1))
      colnames(new_sites) <- c("genomicPos")
      new_sites <- new_sites %>% mutate(id=tx,type="extended")
      tot_new <- rbind(tot_new,new_sites)
    }
    y <- rbind(y,tot_new) %>% select(-type) %>% distinct() %>% as.data.frame()
    rownames(y) <- NULL
    return(y)
  })
  x <- bind_rows(x)
  return(x)
})
peakcalled_1_iso <- get_presence(peakcalled_1_iso) %>%
  left_join(select(assembly,id,ORF) ,by="id") %>%
  select(-id) %>%
  dplyr::relocate(ORF, .after=genomicPos)

peakcalled_multiple_iso <- lapply(file_list,function(x){
  linename <- str_remove(basename(x),paste0("_modified_sites_",which_nucl,".xls"))
  x <- read.xlsx(x,sheetName =paste0("peakcalled_canonical_",which_nucl,"s"),as.data.frame = T)%>%
    mutate(cell_type=linename) %>%
    subset(id %in% assembly_multiple_isoform$id)%>%
    select(-fleming_presence)
  x <- split(x, x$ORF)
  x <- lapply(x, function(y){
    ORF_id <- unique(y$ORF)
    y <-y %>% 
      select("genomicPos","ORF")%>%
      distinct()%>%
      mutate(type="central")
    tot_new <- data.frame()
    for (row in 1:nrow(y)){
      genpos <- y[row, "genomicPos"]
      new_sites <- as.data.frame(seq((genpos-4),(genpos+4),1))
      colnames(new_sites) <- c("genomicPos")
      new_sites <- new_sites %>% mutate(ORF=ORF_id,type="extended")
      tot_new <- rbind(tot_new,new_sites)
    }
    y <- rbind(y,tot_new) %>% select(-type) %>% distinct() %>% as.data.frame()
    rownames(y) <- NULL
    return(y)
  })
  x <- bind_rows(x)
  return(x)
})
peakcalled_multiple_iso <- get_presence(peakcalled_multiple_iso)
  
for (i in 2:3){
  shared <- rbind(peakcalled_1_iso,peakcalled_multiple_iso) %>% 
    subset(presence>=i) %>% 
    distinct()
  shared <- split(shared,shared$ORF)
  shared <- lapply(X=shared, function(x){get_shared(x)})
  shared <- bind_rows(shared) %>%
    rowwise()%>%
    dplyr::mutate(fleming_presence=get_fleming(left_interval,right_interval)) %>%
    rowwise()%>%
    dplyr::mutate(fragments_partial_included=get_fragment_partial_overlap(left_interval,right_interval))%>%
    rowwise()%>%
    dplyr::mutate(nucleotides_in_junctions=get_ORF_per_nucl(ORF,left_interval,right_interval))
  if(i==2){shared_2<-shared}
  else{shared_3 <- shared}
}


write.xlsx(as.data.frame(shared_2), file=rdp(paste0("shared_",which_nucl,".xls")),sheetName=paste0("shared_2_",which_nucl),row.names=F,col.names=T)
write.xlsx(as.data.frame(shared_3), file=rdp(paste0("shared_",which_nucl,".xls")),sheetName=paste0("shared_3_",which_nucl),row.names=F,col.names=T,append=T)


############################### TRACKS #########################################

# Canonical assembly track
assembly_track <- assembly %>%
  subset(canonicity=="C") %>%
  select(-canonicity,-id)%>%
  relocate(ORF, .after = end)%>%
  select(-others)
write.table(assembly_track,sep = "\t",quote = F,row.names = F,col.names = F,file = rdp("shared/bedtracks/assembly.bed"))

# mods

final2 <- subset(shared_2,mod!="")
final2 <- split(final2,final2$ORF)
bed_list <- lapply(final2, function(x){
  ORF_id <- unique(x$ORF)
  x <- x %>%
    mutate(start=(left_interval-1)) %>%
    mutate(end=right_interval) %>%
    mutate(name=paste0("U_genomicpos=",left_interval,"_",right_interval,"_ORF=",ORF_id))%>%
    mutate(chr="NC_045512v2")%>%
    mutate(score=0)%>%
    mutate(strand="+") %>%
    select(chr,start,end,name,score,strand)
  write.table(x,sep = "\t",quote = F,row.names = F,col.names = F,file = rdp(paste("shared/bedtracks/2_shared_peakcalled_Us_",ORF_id,".bed",sep="")))
  
})

final3 <- subset(shared_3,mod!="")
final3 <- split(final3,final3$ORF)
bed_list <- lapply(final3, function(x){
  ORF_id <- unique(x$ORF)
  x <- x %>%
    mutate(start=(left_interval-1)) %>%
    mutate(end=right_interval) %>%
    mutate(name=paste0("U_genomicpos=",left_interval,"_",right_interval,"_ORF=",ORF_id))%>%
    mutate(chr="NC_045512v2")%>%
    mutate(score=0)%>%
    mutate(strand="+") %>%
    select(chr,start,end,name,score,strand)
  write.table(x,sep = "\t",quote = F,row.names = F,col.names = F,file = rdp(paste("shared/bedtracks/3_shared_peakcalled_Us_",ORF_id,".bed",sep="")))
  
})



############################### PLOTS #########################################

##############  Violin plots
### inspect Fleming sites
ref_id_others <- ref_id_others %>% dplyr::rename(id=ref_id)
flem <- shared_2 %>%
  subset(fleming_presence==T)%>%
  subset(ORF %!in% assembly_multiple_isoform$name) %>%
  mutate(id=get_tx(ORF)) %>%
  left_join(ref_id_others,by="id")
tx_id <- paste0(unique(flem$ref_id),"::",unique(flem$others))
tx_sub <- sub("\\|","_",sub(":","_",sub("::","_",tx_id)))

pdf(paste0("/Users/camillaugolini/Desktop/",flem$left_interval,"-",flem$right_interval,".pdf"),height=20,width=10)
for (cell_line in c("caco2","vero","calu3")){
  if (cell_line=="vero" | cell_line=="caco2"){
    list_tx <- read_tsv(bdp(paste0("analysis/per_cell_line/",cell_line,"/NANOCOMPORE/sampcomp/",tx_sub,"/outnanocompore_results.tsv")),col_types = "ncncccnnncncnc") 
  }
  else if (cell_line=="calu3"){
    list_tx <- read_tsv(ddp2(paste0("sraf_calu3/",tx_sub,"/out_nanocompore_results.tsv")),col_types = "ncncccnnncncnc") 
  }
  GMM_plot <- list_tx %>% 
    subset(genomicPos>=flem$left_interval & genomicPos<=flem$right_interval) %>%
    mutate(Color = ifelse((genomicPos %in% fleming_paper_notation)==T, "red", "black")) %>%
    ggplot(.,aes(x=as.character(genomicPos),y=-log10(GMM_logit_pvalue), color=Color)) +
    geom_point()+
    scale_color_identity()+
    ggtitle(paste0("Site ",flem$left_interval,"-",flem$right_interval,", ",cell_line),subtitle = str_remove(tx_id,"::.*"))+
    vlab
  KS_int_plot<-list_tx %>% 
    subset(genomicPos>=flem$left_interval & genomicPos<=flem$right_interval) %>%
    mutate(Color = ifelse((genomicPos %in% fleming_paper_notation)==T, "red", "black")) %>%
    ggplot(.,aes(x=as.character(genomicPos),y=-log10(KS_intensity_pvalue), color=Color)) +
    geom_point()+
    scale_color_identity()+
    ggtitle(paste0("Site ",flem$left_interval,"-",flem$right_interval,", ",cell_line),subtitle = str_remove(tx_id,"::.*"))+
    vlab
  KS_dwell_plot<-list_tx %>%
    subset(genomicPos>=(flem$left_interval-10) & genomicPos<=(flem$right_interval+10)) %>%
    mutate(Color = ifelse((genomicPos %in% fleming_paper_notation)==T, "red", "black")) %>%
    ggplot(.,aes(x=as.character(genomicPos),y=-log10(KS_dwell_pvalue), color=Color)) +
    geom_point()+
    scale_color_identity()+
    ggtitle(paste0("Site ",flem$left_interval,"-",flem$right_interval,", ",cell_line),subtitle = str_remove(tx_id,"::.*"))+
    vlab
  print(ggarrange(GMM_plot,KS_int_plot,KS_dwell_plot,ncol = 1,nrow = 3))
  }
dev.off()

## inspect UGUAR sites

motifs <- c("UGUAG","UGUAA","UCUAA","UAUAA")
site <- shared_3 %>%
  subset(grepl("UGUAG",mod)==T | grepl("UGUAA",mod)==T | grepl("UCUAA",mod)==T | grepl("UAUAA",mod)==T)%>%
  subset(ORF %!in% assembly_multiple_isoform$name) %>%
  mutate(ref_id=get_tx(ORF)) %>%
  left_join(ref_id_others,by="ref_id")
for (row in 1:nrow(site)){
  tx_id <- paste0(unique(site[row,"ref_id"]),"::",unique(site[row,"others"]))
  tx_sub <- sub("\\|","_",sub(":","_",sub("::","_",tx_id)))
  my_var <- as.data.frame(site[row,])
  pdf(paste0("/Users/camillaugolini/Desktop/",my_var$left_interval,"-",my_var$right_interval,".pdf"),height=20,width=14)
  for (cell_line in c("caco2","vero","calu3")){
    if (cell_line=="vero" | cell_line=="caco2"){
      list_tx <- read_tsv(bdp(paste0("analysis/per_cell_line/",cell_line,"/NANOCOMPORE/sampcomp/",tx_sub,"/outnanocompore_results.tsv")),col_types = "ncncccnnncncnc") 
    }
    else if (cell_line=="calu3"){
      list_tx <- read_tsv(ddp2(paste0("sraf_calu3/",tx_sub,"/out_nanocompore_results.tsv")),col_types = "ncncccnnncncnc") 
    }
    GMM_plot <- list_tx %>% 
      subset(genomicPos>=(my_var$left_interval-1-5) & genomicPos<=(my_var$right_interval-1+5)) %>%
      mutate(x_axis=paste0(genomicPos,"\n",ref_kmer)) %>%
      mutate(Color = ifelse((gsub("T","U",ref_kmer) %in% motifs)==T, "red", "black")) %>%
      ggplot(.,aes(x=as.character(x_axis),y=-log10(GMM_logit_pvalue), color=Color)) +
      geom_point()+
      scale_color_identity()+
      ggtitle(paste0("Site ",my_var$left_interval,"-",my_var$right_interval,", ",cell_line),subtitle = str_remove(tx_id,"::.*"))
    KS_int_plot <- list_tx %>% 
      subset(genomicPos>=(my_var$left_interval-1-5) & genomicPos<=(my_var$right_interval-1+5)) %>%
      mutate(x_axis=paste0(genomicPos,"\n",ref_kmer)) %>%
      mutate(Color = ifelse((gsub("T","U",ref_kmer) %in% motifs)==T, "red", "black")) %>%
      ggplot(.,aes(x=as.character(x_axis),y=-log10(KS_intensity_pvalue), color=Color)) +
      geom_point()+
      scale_color_identity()+
      ggtitle(paste0("Site ",my_var$left_interval,"-",my_var$right_interval,", ",cell_line),subtitle = str_remove(tx_id,"::.*"))
    KS_dwell_plot <- list_tx %>% 
      subset(genomicPos>=(my_var$left_interval-1-10) & genomicPos<=(my_var$right_interval-1+10)) %>%
      mutate(x_axis=paste0(genomicPos,"\n",ref_kmer)) %>%
      mutate(Color = ifelse((gsub("T","U",ref_kmer) %in% motifs)==T, "red", "black")) %>%
      ggplot(.,aes(x=as.character(x_axis),y=-log10(KS_dwell_pvalue), color=Color)) +
      geom_point()+
      scale_color_identity()+
      ggtitle(paste0("Site ",my_var$left_interval,"-",my_var$right_interval,", ",cell_line),subtitle = str_remove(tx_id,"::.*"))+
      vlab
    print(ggarrange(GMM_plot,KS_int_plot,KS_dwell_plot,ncol = 1,nrow = 3))
  }
  dev.off()
}



 

ass <- assembly[3,] %>% 
  separate(blockSizes, sep=",", into=c("block_size_1","block_size_2")) %>%
  separate(blockStarts, sep=",", into=c("block_start_1","block_start_2"))
start_block_1 <- ass$start
end_block_1 <- as.numeric(ass$start)+as.numeric(ass$block_size_1)
genomic_block_start_2 <- as.numeric(ass$start)+as.numeric(ass$block_start_2)
if(my_cord<=end_block_1 & my_cord>=start_block_1) {
  trascr_cord <- my_cord-start_block_1-1
}else if(my_cord>=end_block_1){
  trascr_cord <- my_cord-1-genomic_block_start_2+as.numeric(ass$block_size_1)
}




x <- final[[9]]
full_seq <- seq(min(x$genomicPos),max(x$genomicPos),1) %>%
  as.data.frame()
colnames(full_seq)<- c("genomicPos")
x<-as.data.frame(left_join(full_seq,x,by="genomicPos"))
right <- data.frame()
if(is.na(x[2,"ORF"])==T){right <- as.data.frame(x[1,"genomicPos"])}
left <- as.data.frame(x[1,"genomicPos"])
for (row in 2:nrow(x)){
  if(is.na(x[(row+1),"ORF"])==T & is.na(x[(row-1),"ORF"])==F & is.na(x[(row),"ORF"])==F){
    right <- rbind(right,x[row,"genomicPos"])}
  else if(is.na(x[(row+1),"ORF"])==T & is.na(x[(row-1),"ORF"])==T & is.na(x[(row),"ORF"])==F){
    right <- rbind(right,x[row,"genomicPos"])}
}
colnames(right) <- c("right_interval")
colnames(left) <- c("left_interval")
for (row in 2:nrow(x)){
  if(is.na(x[(row+1),"ORF"])==F & is.na(x[(row-1),"ORF"])==T & is.na(x[(row),"ORF"])==F){
    left <- rbind(left,x[row,"genomicPos"])}
  else if(is.na(x[(row+1),"ORF"])==T & is.na(x[(row-1),"ORF"])==T & is.na(x[(row),"ORF"])==F){
    left <- rbind(left,x[row,"genomicPos"])}
}
colnames(left) <- c("left_interval")
interval <- cbind(left,right)














sites_all <- lapply(file_list,function(x){
  linename <- str_remove(basename(x),"_modified_sites.xls")
  x <- read.xlsx(x,sheetName ="all_significant_sites",as.data.frame = T)
  x <- x %>%
    select(pos,chr,genomicPos,ref_id,ref_kmer,ORF,canonicity,motif,junction_type,fragment_ID,fleming_presence)%>%
    mutate(cell_type=linename)
})
cols_selection <- c("pos","chr","genomicPos","ref_id","ref_kmer","ORF","canonicity","motif","junction_type","fragment_ID","fleming_presence")

sites_all_1_isoform <- lapply(sites_all, function(x){x %>% select(-cell_type)})
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
  select(chr,genomicPos,ref_kmer,ORF,canonicity,motif,junction_type,fragment_ID,fleming_presence)%>%
  mutate(pos="NA")%>%
  mutate(ref_id="NA")

sites_all<- rbind(sites_all_1_isoform,sites_all_multiple_isoform)               # all significant sites                                
canonical_signif_U <- sites_all %>%                                             # canonical U sites (no junction)
  subset(junction_type=="No junction") %>% 
  subset(substring(ref_kmer,3,3)=="U")%>%
  subset(canonicity=="C")
signif_5p <- sites_all %>%                                                      # significant sites before 100 
  subset((genomicPos-3)<=100)
signif_5p_canonical_U <- sites_all %>%                                          # significant canonical U sites before 100 
  subset((genomicPos-3)<=100) %>%
  subset(canonicity=="C") %>%
  subset(substring(ref_kmer,3,3)=="U")



file_list <- list.files(path = RESULTSDIR,pattern = "*_modified_sites.xls" , full.names = TRUE,  recursive = F)
peakcalled <- lapply(file_list,function(x){
  linename <- str_remove(basename(x),"_modified_sites.xls")
  x <- read.xlsx(x,sheetName ="peakcalled_canonical_Us",as.data.frame = T)%>%
    mutate(cell_type=linename) %>%
    subset(ref_id %!in% assembly_multiple_isoform$id)
  x <- split(x, x$ref_id)
  x <- lapply(x, function(y){
    WD=rdp(paste("shared/tracks_for_overlap/",gsub("\\|","_", unique(y$ref_id)),sep=""))
    dir.create(WD, showWarnings = FALSE)
    y <-y %>% 
      select("chr", "genomicPos","ref_id")%>%
      mutate(start=(genomicPos-6),end=(genomicPos+6),id=paste(ref_id,start,end,sep='_'))%>%
      select(-ref_id,-genomicPos) %>% 
      mutate(score="*",strand="+")
    write.table(y,file = paste0(WD,"/",linename,".bed"), sep = "\t",quote = F,row.names = F,col.names = F)
    return(y)
  })
})

dir_list <- list.dirs(path = rdp("shared/tracks_for_overlap"), full.names = TRUE,  recursive = F)
shared_3 <- lapply(dir_list,function(x){
  tx <- basename(x)
  x <- read.table(file=paste0(x,"/common_merged.bed3"),sep = '\t')
  colnames(x) <- c("chr","start","end")

})







cols_selection <- c("pos","chr","genomicPos","ref_id","strand","ORF","canonicity","junction_type","fragment_ID","fleming_presence")

peakcalled_1_isoform <- peakcalled %>%
  bind_rows() %>%
  subset(ref_id %!in% assembly_multiple_isoform$id) %>%
  distinct()
peakcalled_1_isoform <- split(peakcalled_1_isoform,peakcalled_1_isoform$ref_id)
bed_list <- lapply(peakcalled_1_isoform, function(x){
  tx_id <- unique(x$ref_id)
  x <- x %>%
    mutate(start=(genomicPos-6)) %>%
    mutate(end=(genomicPos+6)) %>%
    mutate(name=paste0("genomicpos=",genomicPos,"_tx=",tx_id))%>%
    mutate(score=0)%>%
    mutate(strand="+") %>%
    select(chr,start,end,name,score,strand)
  write.table(x,sep = "\t",quote = F,row.names = F,col.names = F,file = rdp(paste("shared/tracks_for_overlap/",cell_line,"_peakcalled_Us_",tx_id,".bed",sep="")))
  
})



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
  dplyr::mutate(presence=sum(c(present_in_caco2,present_in_calu3,present_in_vero), na.rm = TRUE))%>%
  as.data.frame()

# Fleming

fleming_redundant_signif <- sites_all %>%                                       # all Fleming significant sites
  subset(fleming_presence==T) 
fleming_nonredundant_signif <- sites_all %>%                                  # all Fleming significant non redundant sites
  subset(fleming_presence==T) %>%
  select(genomicPos,ref_kmer)%>%
  dplyr::group_by(genomicPos,ref_kmer) %>% 
  dplyr::filter(row_number() == 1)


fleming_all <- lapply(file_list,function(x){
  linename <- str_remove(basename(x),"_modified_sites.xls")
  x <- read.xlsx(x,sheetName ="Fleming_redundant_ALL_sites",as.data.frame = T)
  x <- x %>%
    select(pos,chr,genomicPos,ref_id,ref_kmer,ORF,canonicity,motif,junction_type,fragment_ID,fleming_presence)%>%
    mutate(cell_type=linename)
})
cols_selection <- c("pos","chr","genomicPos","ref_id","ref_kmer","ORF","canonicity","motif","junction_type","fragment_ID","fleming_presence")

fleming_all_1_isoform <- lapply(fleming_all, function(x){x%>%select(-cell_type)})
fleming_all_1_isoform <-plyr::join_all(fleming_all_1_isoform, by=cols_selection,type="inner") %>%
  subset(ref_id %!in% assembly_multiple_isoform$id)

fleming_all_multiple_isoform <- lapply(fleming_all, function(x){
  x<- x %>%
    subset(ref_id %in% assembly_multiple_isoform$id) 
})
fleming_all_multiple_isoform <- bind_rows(fleming_all_multiple_isoform)
fleming_all_multiple_isoform <- split(fleming_all_multiple_isoform,fleming_all_multiple_isoform$genomicPos)
fleming_all_multiple_isoform <- lapply(fleming_all_multiple_isoform,function(x){
  x <-split(x,x$ORF)
  x <- lapply(x,function(y){
    if(length(unique(y$cell_type))==3){y<-y%>%slice(1:1)%>%mutate(del=F)}
    else{y<-y%>%mutate(del=T)}
  })
  x <- bind_rows(x)
})
fleming_all_multiple_isoform <-bind_rows(fleming_all_multiple_isoform) %>%
  subset(del==F)%>%
  select(chr,genomicPos,ref_kmer,ORF,canonicity,motif,junction_type,fragment_ID,fleming_presence)%>%
  mutate(pos="NA")%>%
  mutate(ref_id="NA")

fleming_all<- rbind(fleming_all_1_isoform,fleming_all_multiple_isoform)         # all Fleming sites
fleming_nonredundant_all <- fleming_all %>%                                     # all Fleming non redundant sites
  dplyr::select(genomicPos,ref_kmer)%>%
  dplyr::group_by(genomicPos,ref_kmer) %>% 
  dplyr::filter(row_number() == 1)
fleming_nonredundant_all <- as.data.frame(fleming_nonredundant_all)





############################# WRITE TABLES #####################################

# Table with identity of shared modifications

write.xlsx(sites_all, rdp("shared_sites.xls"),sheetName="all_significant_sites",row.names=F,col.names=T)
write.xlsx(canonical_signif_U, rdp("shared_sites.xls"),sheetName="all_canonical_Us_sign_sites",row.names=F,col.names=T,append=TRUE)
write.xlsx(signif_5p, rdp("shared_sites.xls"),sheetName="5p_significant_sites",row.names=F,col.names=T,append=TRUE)
write.xlsx(signif_5p_canonical_U, rdp("shared_sites.xls"),sheetName="5p_canonical_Us_sign_sites",row.names=F,col.names=T,append=TRUE)
write.xlsx(final, rdp("shared_sites.xls"),sheetName="peakcalled_canonical_Us",row.names=F,col.names=T,append=TRUE)

# Table with Fleming sites
write.xlsx(fleming_all, rdp("shared_sites.xls"),sheetName="Fleming_redundant_ALL_sites",row.names=F,col.names=T,append=TRUE)
write.xlsx(fleming_redundant_signif, rdp("shared_sites.xls"),sheetName="Fleming_redundant_sign_sites",row.names=F,col.names=T,append=TRUE)
write.xlsx(fleming_nonredundant_signif, rdp("shared_sites.xls"),sheetName="Fleming_non_red_CandNC_sign_sites",row.names=F,col.names=T,append=TRUE)
write.xlsx(fleming_nonredundant_all, rdp("shared_sites.xls"),sheetName="Fleming_non_red_CandNC_all_sites",row.names=F,col.names=T,append=TRUE)



###  Working on high-confidence pseudoU sites

f <- list.files(rdp("shared/corresp_tables"), include.dirs = F, full.names = T, recursive = T)# remove the files
file.remove(f)

extract_others <- tx %>%
  select(name) %>%
  separate(name,into=c("ref_id","others"),sep="::") %>%
  separate(others,into=c("others","del"),sep="#") %>%
  select(-del)
extract_kmers <- sites_all_1_isoform %>%
  select(pos,genomicPos,ref_id,ref_kmer)
tables_1_iso <- peakcalled_1_isoform %>%
  as.data.frame()%>%
  rowwise()%>%
  mutate(presence=sum(c(present_in_caco2,present_in_calu3,present_in_vero), na.rm = TRUE))%>%
  as.data.frame() %>%
  subset(presence==3) %>%
  left_join(extract_kmers,by=c("ref_id","pos","genomicPos")) %>%
  mutate(genomicPos=(genomicPos-3)) %>% 
  left_join(extract_others,by="ref_id") %>%
  mutate(ref_id=paste0(ref_id,"::",others))

tables_1_iso <- split(tables_1_iso,tables_1_iso$genomicPos)
lapply(tables_1_iso,function(x){
  genPos <- unique(x$genomicPos)
  ref_kmer <- unique(x$ref_kmer)
  x <- x %>%
    select(pos,ref_id)
  write.table(x,sep = "\t",quote = F,row.names = F,col.names = F,file = rdp(paste("shared/corresp_tables/",ref_kmer,"_",genPos,".txt",sep="")))
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




# Working on high-confidence sites

corresp_table <- read.xlsx(rdp(paste0(cell_line,"_modified_sites.xls")),sheetName="Fleming_redundant_ALL_sites",as.data.frame = T)
corresp_table <- split(corresp_table,corresp_table$genomicPos)
lapply(corresp_table,function(x){
  genPos <- unique(x$genomicPos)
  ref_kmer <- unique(x$ref_kmer)
  x <- x %>%
    mutate(ref_tx=paste(ref_id,"::",others,sep = ""))%>%
    select(pos,ref_tx)
  write.table(x,sep = "\t",quote = F,row.names = F,col.names = F,file = paste(CORRESPTABLEDIR,"/",ref_kmer,"_",genPos,".txt",sep=""))
})

