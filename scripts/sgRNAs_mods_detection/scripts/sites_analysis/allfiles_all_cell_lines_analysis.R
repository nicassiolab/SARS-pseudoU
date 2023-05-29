library(rstudioapi)

source(paste(dirname(getSourceEditorContext()$path),"variables.R",sep="/"))
source(paste(dirname(getSourceEditorContext()$path),"functions.R",sep="/"))
source(paste(dirname(getSourceEditorContext()$path),"additional_dataframes.R",sep="/"))

dir.create(paste0(RESULTSDIR,"/shared/bedtracks"),showWarnings = F,recursive = T)
dir.create(paste0(RESULTSDIR,"/shared/corresp_tables"),showWarnings = F,recursive=T)
dir.create(paste0(RESULTSDIR,"/shared/tracks_for_overlap"),showWarnings = F,recursive=T)


####################### PROCESSING OF THE DATABASES ############################

file_list <- list.files(path = RESULTSDIR,
                        pattern = paste0("*_modified_sites_",which_nucl,".xls"), 
                        full.names = TRUE,
                        recursive = F)


# For every significant k-mer, detect its presence in each cell line  peakcalled dataset

# start with sgRNAs with a unique isoform
peakcalled_1_iso <- lapply(file_list,function(x){
  linename <- str_remove(basename(x),paste0("_modified_sites_",which_nucl,".xls"))
  x <- read.xlsx(x,sheetName = paste0("peakcalled_canonical_",which_nucl,"s"),
                 as.data.frame = T) %>%
    mutate(cell_type=linename) %>%
    subset(id %!in% assembly_multiple_isoform$id) %>%
    dplyr::select(-fleming_presence)
  x <- split(x, x$id)
  x <- lapply(x, function(y){
    tx <- unique(y$id)
    y <-y %>%
      dplyr::select("genomicStartKmer","id") %>%
      mutate(genomicPos=(genomicStartKmer+2)) %>%
      mutate(type="central") %>%
      dplyr::select(-genomicStartKmer)
    tot_new <- data.frame()
    for (row in 1:nrow(y)){
      genpos <- y[row, "genomicPos"]
      new_sites <- as.data.frame(seq((genpos-4),(genpos+4),1))
      colnames(new_sites) <- c("genomicPos")
      new_sites <- new_sites %>% mutate(id=tx,type="extended")
      tot_new <- rbind(tot_new,new_sites)
    }
    y <- rbind(y,tot_new) %>%
      dplyr::select(-type) %>%
      distinct() %>%
      as.data.frame()
    rownames(y) <- NULL
    return(y)
  })
  x <- bind_rows(x)
  return(x)
})
peakcalled_1_iso <- get_presence(peakcalled_1_iso) %>%
  left_join(dplyr::select(names,id,name), by="id") %>%
  dplyr::select(-id) %>%
  dplyr::relocate(name, .after=genomicPos) %>%
  dplyr::rename(ORF=name)


# then process sgRNAs with multiple isoforms
peakcalled_multiple_iso <- lapply(file_list,function(x){
  linename <- str_remove(basename(x),paste0("_modified_sites_",which_nucl,".xls"))
  x <- read.xlsx(x,sheetName =paste0("peakcalled_canonical_",which_nucl,"s"),
                 as.data.frame = T)%>%
    mutate(cell_type=linename) %>%
    subset(id %in% assembly_multiple_isoform$id) %>%
    dplyr::select(-fleming_presence) 
  x <- split(x, x$ORF)
  x <- lapply(x, function(y){
    ORF_id <- unique(y$ORF)
    y <-y %>%
      mutate(genomicPos=(genomicStartKmer+2)) %>%
      dplyr::select(-genomicStartKmer) %>%
      dplyr::select("genomicPos","ORF") %>%
      distinct() %>%
      mutate(type="central")
    tot_new <- data.frame()
    for (row in 1:nrow(y)){
      genpos <- y[row, "genomicPos"]
      new_sites <- as.data.frame(seq((genpos-4),(genpos+4),1))
      colnames(new_sites) <- c("genomicPos")
      new_sites <- new_sites %>%
        mutate(ORF=ORF_id,type="extended")
      tot_new <- rbind(tot_new,new_sites)
    }
    y <- rbind(y,tot_new) %>%
      dplyr::select(-type) %>%
      distinct() %>%
      as.data.frame()
    rownames(y) <- NULL
    return(y)
  })
  x <- bind_rows(x)
  return(x)
})
peakcalled_multiple_iso <- get_presence(peakcalled_multiple_iso)
  

# Join single-isoform dataset and multiple-isoform dataset and get sites present in at least 2 or 3 cell lines
for (i in 2:3){
  shared <- rbind(peakcalled_1_iso,peakcalled_multiple_iso) %>% 
    subset(presence>=i) %>% 
    distinct()
  shared <- split(shared,shared$ORF)
  shared <- lapply(X=shared, function(x){get_shared(x)})
  shared <- bind_rows(shared) %>%
    rowwise()%>%
    dplyr::mutate(fleming_presence=get_fleming(left_interval,right_interval)) %>%
    rowwise() %>%
    dplyr::mutate(fragments_partial_included=get_fragment_partial_overlap(left_interval,right_interval)) %>%
    rowwise() %>%
    dplyr::mutate(nucleotides_in_junctions=get_ORF_per_nucl(ORF,left_interval,right_interval))
  if(i==2){shared_2<-shared}
  else{shared_3 <- shared}
}


shared_2_bed <- shared_2 %>%
  dplyr::rename(bed_start=left_interval,
                bed_end=right_interval,
                nucleotides_in_junctions_0_based=nucleotides_in_junctions) %>%
  mutate(bed_end=(bed_end+1))

shared_3_bed <- shared_3 %>%
  dplyr::rename(bed_start=left_interval,
                bed_end=right_interval,
                nucleotides_in_junctions_0_based=nucleotides_in_junctions) %>%
  mutate(bed_end=(bed_end+1))

write.xlsx(as.data.frame(shared_2_bed), file=rdp(paste0("shared_",which_nucl,".xls")),
           sheetName=paste0("shared_2_",which_nucl),row.names=F,col.names=T)
write.xlsx(as.data.frame(shared_3_bed), file=rdp(paste0("shared_",which_nucl,".xls")),
           sheetName=paste0("shared_3_",which_nucl),row.names=F,col.names=T,append=T)



############################### TRACKS #########################################

# Part of the code to generate tracks for the UCSC Genome Browser

# Canonical assembly track
assembly_track <- assembly_named_canonical %>%
  dplyr::rename(ORF=name) %>%
  mutate(id=paste0(id,"_",ORF)) %>%
  dplyr::select(-canonicity,-others,-uniprot_id,-tip,-qip)%>%
  relocate(id, .after = end)
write.table(assembly_track,sep = "\t",quote = F,row.names = F,col.names = F,file = rdp("shared/bedtracks/assembly.bed"))

# Track with sites present in at least two cell lines

final2 <- as.data.frame(shared_2_bed) %>%
  subset(mod!="") %>%
  rowwise() %>% 
  mutate(UNUAR_presence=get_UNUAR(mod)) %>%
  subset(UNUAR_presence==T) %>%
  dplyr::select(-UNUAR_presence)


final2 <- split(final2,final2$ORF)
bed_list <- lapply(final2, function(x){
  ORF_id <- unique(x$ORF)
  x <- x %>%
    dplyr::rename(start=bed_start,end=bed_end) %>%
    mutate(name=paste0("U_genomicpos=",start,"_",end,"_ORF=",ORF_id))%>%
    mutate(chr="NC_045512v2")%>%
    mutate(score=0)%>%
    mutate(strand="+") %>%
    dplyr::select(chr,start,end,name,score,strand)
  write.table(x,sep = "\t",quote = F,row.names = F,col.names = F,file = rdp(paste("shared/bedtracks/2_shared_peakcalled_Us_",ORF_id,".bed",sep="")))
  
})

# Track with sites present in three cell lines

final3 <- as.data.frame(shared_3_bed) %>%
  subset(mod!="") %>%
  rowwise() %>% 
  mutate(UNUAR_presence=get_UNUAR(mod)) %>%
  subset(UNUAR_presence==T) %>%
  dplyr::select(-UNUAR_presence)
final3 <- split(final3,final3$ORF)
bed_list <- lapply(final3, function(x){
  ORF_id <- unique(x$ORF)
  x <- x %>%
    dplyr::rename(start=bed_start,end=bed_end) %>%
    mutate(name=paste0("U_genomicpos=",start,"_",end,"_ORF=",ORF_id))%>%
    mutate(chr="NC_045512v2")%>%
    mutate(score=0)%>%
    mutate(strand="+") %>%
    dplyr::select(chr,start,end,name,score,strand)
  write.table(x,sep = "\t",quote = F,row.names = F,col.names = F,file = rdp(paste("shared/bedtracks/3_shared_peakcalled_Us_",ORF_id,".bed",sep="")))
  
})


