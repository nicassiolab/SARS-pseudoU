library(tidyverse)
library(tidyr)
library(R.utils)
library(reshape2)
library(rlang)
library(stringr)
library(plyr)
library(dplyr)
library(xlsx)
library(seqinr)
library(ggpubr)



ROOTDIR="/Volumes/scratch/TSSM/cugolini/cov"
RESULTSDIR="/Volumes/scratch/FN/TL/cugolini/cov/scripts/downstream/results_allfiles_LOR05_pval001"
DATADIR="/Volumes/scratch/FN/TL/cugolini/cov/analysis"
DATADIR2="/Volumes/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/nanocompore/comparison"
dir.create(paste0(RESULTSDIR,"/shared/bedtracks"),showWarnings = F)
dir.create(paste0(RESULTSDIR,"/shared/corresp_tables"),showWarnings = F)
dir.create(paste0(RESULTSDIR,"/shared/tracks_for_overlap"),showWarnings = F)


########################### PARAMETERS ##########################################

`%!in%` <- Negate(`%in%`)
fleming_paper_notation <- c(22322,23317,27164,28417,28759,28927,29418)


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
get_motif = function(seq_to_search) {
  motifs <- sapply(as.vector(sitelist$ref_kmer), function(x){
    x <-str_extract_all(seq_to_search,x)
  })
  motifs <-toString(unlist(motifs,use.names = F))
  return(motifs)
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



########################### GENERAL DATA #######################################

tx <- read_tsv(bdp("analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/orf_annotate/orf_annotate.bed"), col_names=c("chr", "start", "end", "name", "score", "strand", "cdsStart", "cdsEnd", ".", "ex", "exLen", "exSt"), col_types="cnncncnnnncc")
ref_id_others <- tx %>% 
  select(name) %>% 
  separate(name, into=c("ref_id","others"),sep="::") %>% 
  separate(others, into=c("others","bla"),sep="#") %>% 
  select(-bla)
viral_genome <- read.fasta(file="/Volumes/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa") %>%
  getSequence(as.string=T) %>% unlist()
assembly <-read.table(bdp("scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus.bed"), col.names = c("chrom","start","end","id","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"), sep="\t")
sitelist <- read_tsv(bdp("scripts_new/backupped_data/sites_pseudoU_motifs_blacklist.txt"),col_types="cccc") %>% 
  dplyr::rename(motif=IUPAC) 
file_list <- list.files(path = RESULTSDIR,pattern = "*_modified_sites.xls" , full.names = TRUE,  recursive = F)


# Dataframe that reports canonicity for every transcript of the assembly
canonicity <- dplyr::select(assembly,start,end,id) %>%
  mutate(canonicity=ifelse(start>100,"NC", "C")) %>%
  mutate(canonicity=ifelse(end<29000,"NC", canonicity)) %>%
  #mutate(canonicity=ifelse(id=="efad7b96-ac2e-4ce1-9b83-863ffdb18eac|116","NC", canonicity)) %>%  #ORF10
  #mutate(canonicity=ifelse(id=="de81ef19-655d-4ced-a9cc-cb8384001058|107","NC", canonicity)) %>%  #ORF9D
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

### peak called datasets

peakcalled_1_iso <- lapply(file_list,function(x){
  linename <- str_remove(basename(x),"_modified_sites.xls")
  x <- read.xlsx(x,sheetName ="peakcalled_canonical_Us",as.data.frame = T)%>%
    mutate(cell_type=linename) %>%
    subset(ref_id %!in% assembly_multiple_isoform$id)
  x <- split(x, x$ref_id)
  x <- lapply(x, function(y){
    tx <- unique(y$ref_id)
    y <-y %>% 
      select("genomicPos","ref_id")%>%
      mutate(type="central")
    tot_new <- data.frame()
    for (row in 1:nrow(y)){
      genpos <- y[row, "genomicPos"]
      new_sites <- as.data.frame(seq((genpos-4),(genpos+4),1))
      colnames(new_sites) <- c("genomicPos")
      new_sites <- new_sites %>% mutate(ref_id=tx,type="extended")
      tot_new <- rbind(tot_new,new_sites)
    }
    y <- rbind(y,tot_new) %>% select(-type) %>% distinct()
    return(y)
  })
  x <- bind_rows(x)
})

appo <- peakcalled_1_iso
peakcalled_1_iso <- bind_rows(peakcalled_1_iso)
present_in_caco2 <- do.call(paste0, peakcalled_1_iso) %in% do.call(paste0, appo[[1]])
present_in_calu3 <- do.call(paste0, peakcalled_1_iso) %in% do.call(paste0, appo[[2]])
present_in_vero <- do.call(paste0, peakcalled_1_iso) %in% do.call(paste0, appo[[3]])
peakcalled_1_iso <- cbind(peakcalled_1_iso,present_in_caco2,present_in_calu3,present_in_vero)

temp <- assembly %>% select(id,name) %>% dplyr::rename(ref_id=id,ORF=name)
peakcalled_1_iso <- peakcalled_1_iso %>% 
  rowwise() %>%
  dplyr::mutate(presence=sum(c(present_in_caco2,present_in_calu3,present_in_vero), na.rm = TRUE))%>%
  left_join(temp,by="ref_id") %>%
  select(-ref_id) %>%
  dplyr::relocate(ORF, .after=genomicPos)

peakcalled_multiple_iso <- lapply(file_list,function(x){
  linename <- str_remove(basename(x),"_modified_sites.xls")
  x <- read.xlsx(x,sheetName ="peakcalled_canonical_Us",as.data.frame = T)%>%
    mutate(cell_type=linename) %>%
    subset(ref_id %in% assembly_multiple_isoform$id)
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
    y <- rbind(y,tot_new) %>% select(-type) %>% distinct()
    return(y)
  })
  x <- bind_rows(x)
})

appo <- peakcalled_multiple_iso
peakcalled_multiple_iso <- bind_rows(peakcalled_multiple_iso)
present_in_caco2 <- do.call(paste0, peakcalled_multiple_iso) %in% do.call(paste0, appo[[1]])
present_in_calu3 <- do.call(paste0, peakcalled_multiple_iso) %in% do.call(paste0, appo[[2]])
present_in_vero <- do.call(paste0, peakcalled_multiple_iso) %in% do.call(paste0, appo[[3]])
peakcalled_multiple_iso <- cbind(peakcalled_multiple_iso,present_in_caco2,present_in_calu3,present_in_vero)

peakcalled_multiple_iso <- peakcalled_multiple_iso %>% 
  rowwise() %>%
  dplyr::mutate(presence=sum(c(present_in_caco2,present_in_calu3,present_in_vero), na.rm = TRUE))

for (i in 2:3){
  shared <- rbind(peakcalled_1_iso,peakcalled_multiple_iso) %>% 
    subset(presence==i) %>% 
    distinct()
  shared <- split(shared,shared$ORF)
  shared <- lapply(shared, function(x){
    ORF_id<-unique(x$ORF)
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
    interval <- cbind(left,right) %>%
      dplyr::mutate(left_interval=(left_interval-2),right_interval=(right_interval+2))%>%
      rowwise()%>%
      dplyr::mutate(sequence=gsub("T","U",toupper(substr(viral_genome,left_interval, right_interval)))) %>%
      rowwise()%>%
      dplyr::mutate(mod=get_motif(sequence)) %>%
      rowwise()%>%
      dplyr::mutate(fleming_presence=get_fleming(left_interval,right_interval)) %>%
      dplyr::mutate(ORF=ORF_id)
    return(interval)
  })
  shared <- bind_rows(shared)
  if(i==2){
    shared_2<-shared
    }
  else{
    shared_3 <- shared
    }
}

write.xlsx(as.data.frame(shared_2), file="/Users/camillaugolini/Desktop/shared.xls",sheetName="shared_2",row.names=F,col.names=T)
write.xlsx(as.data.frame(shared_3), file="/Users/camillaugolini/Desktop/shared.xls",sheetName="shared_3",row.names=F,col.names=T,append=T)

### inspect Fleming sites

flem <- shared_3 %>%
  subset(fleming_presence==T)%>%
  subset(ORF %!in% assembly_multiple_isoform$name) %>%
  mutate(ref_id=get_tx(ORF)) %>%
  left_join(ref_id_others,by="ref_id")
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
  x%>%
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


# Figure for transcript

ciao <- final %>%
  subset(presence>=2)
ciao <- as.data.frame(table(ciao$genomicPos))
ciao$genomicPos<- as.character(ciao$genomicPos)
ciao$presence<- as.character(ciao$presence)

ggplot(ciao) +
  geom_bar(aes(x = genomicPos, fill = presence))+
  scale_fill_manual(values = c("gray", "red"))+
  vlab

pdf(rdp("FIGURES/supplementary.pdf"), width=10, height=5)
ggplot(data=ciao, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") +
  ylab("N° of transcripts (grouped per ORF) with a modified uridine") +
  xlab("Genomic Position") +
  theme_bw() +
  vlab
dev.off()
  


pippo <- final %>%
  subset(presence>2)%>%
  select(genomicPos)%>%distinct()

pluto <- bind_rows(sites_all) %>%
  select(genomicPos,ref_kmer)%>%distinct()

miao <- left_join(pippo,pluto,by="genomicPos")
miaone <- left_join(miao, peakcalled[[2]])

