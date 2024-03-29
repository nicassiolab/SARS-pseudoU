library(rstudioapi)

source(paste(dirname(getSourceEditorContext()$path),"variables.R",sep="/"))
source(paste(dirname(getSourceEditorContext()$path),"functions.R",sep="/"))
source(paste(dirname(getSourceEditorContext()$path),"additional_dataframes.R",sep="/"))


########################### DIRECTORIES ########################################

dir.create(paste0(RESULTSDIR,"/FIGURES/"),showWarnings = F,recursive = T)
dir.create(paste0(RESULTSDIR,"/BEDTRACKS/"),showWarnings = F,recursive=T)

########################### GENERAL DATA #######################################

# DRS databases available
all_lines_tx <- list.files(path = all_cell_lines_results ,pattern = "*_results.tsv" , 
                           full.names = TRUE,  recursive = T)

####################### PROCESSING OF THE DATABASES ############################


# Extraction of the datasets from Nanocompore data
tx_data <- lapply(all_lines_tx, function(x) {
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
      mutate(significant=ifelse(
        GMM_logit_pvalue<=pval_thresh & abs(Logit_LOR)>=LOR_thresh,T,F)) %>%      # indicate for every sample using TRUE or FALSE if its LOR and pvalue are significant according to the established thresholds
      rowwise()%>%
      mutate(fleming_presence=get_fleming_single_site(genomicPos)) %>%            # add a column to indicate if the site is present in the Fleming et al. paper
      separate(ref_id, into=c("id","others"), sep="::") %>%
      rowwise() %>%
      mutate(ORF = get_ORF(id)) %>%                                               # get the ORF corresponding to the transcript model id
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
    x$fragment_ID <- factor(x$fragment_ID,levels=c("No_Fragment",
                                                   "Fragment1",
                                                   "Fragment1_2",
                                                   "Fragment2_3",
                                                   "Fragment3_4",
                                                   "Fragment4_5",
                                                   "Fragment5",
                                                   "Fragment6",
                                                   "Fragment6_7",
                                                   "Fragment7_8",
                                                   "Fragment8_9",
                                                   "Fragment9_10",
                                                   "Fragment10"))
  return(x)
})  



########################### PLOTS ##############################################


fig_plots <- bind_rows(toplot) %>% subset(canonicity=="C") 
fig_plots <- split(fig_plots,fig_plots$id)
plots <- lapply(fig_plots,function(x){
    p <- x %>%
      {
        ggplot(., aes(x=abs(as.numeric(Logit_LOR)), y=-log10(GMM_logit_pvalue), color=fragment_ID)) +
          geom_point() +
          {if (nrow(subset(x,junction=="No junction" & significant==T & motif=="UNUAR"))>0) 
            ggrepel::geom_label_repel(data=dplyr::filter(.,(junction=="No junction" & significant==T & motif=="UNUAR")),
                                      aes(label=paste0(ref_kmer, " (",genomicPos,")")), 
                                      colour="black", size=7,box.padding = 5, 
                                      min.segment.length = 0,force=2,max.overlaps = Inf)}+
          scale_color_manual(breaks = c("No_Fragment","Fragment1","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10"),
                             values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
          ggtitle(unique(x$ORF), subtitle = paste0(unique(x$id)," ",unique(x$canonicity))) +
          theme_bw()+
          xlab("abs(Log Odds Ratio)")+
          theme(plot.title = element_text(size = 25, face = "bold"), 
                plot.subtitle = element_text(size=25),
                axis.title = element_text(size=25),
                axis.text = element_text(size=15),
                legend.title=element_text(size=20), 
                legend.text=element_text(size=18))+
          geom_hline(yintercept=2, linetype='dotted', col = 'red')+
          geom_vline(xintercept=0.5, linetype='dotted', col = 'red')+ 
          geom_text(aes(0,2,label = "pval_thresh",colour='red', vjust = -1),)
      }
    new_lab<- as.integer((ggplot_build(p)$layout$panel_params[[1]]$y$get_breaks()[2]-ggplot_build(p)$layout$panel_params[[1]]$y$get_breaks()[1])/4)
    p <- p + geom_text(aes(0.5,y=-new_lab,label = "LOR_thresh",colour='red', vjust = -1))
    return(p)
})

pdf(rdp(paste0("/FIGURES/cell_lines_pulled_sharkfin_canonical.pdf")),height=30,width=25)
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
  dplyr::rename(genomicStartKmer = genomicPos) %>%
  mutate(genomicEndKmer = (genomicStartKmer+5)) %>%
  relocate(genomicEndKmer, .after = genomicStartKmer)                            # ALL MODS : significant,no junctions
write.xlsx(final_id_all, rdp(paste0("all_cell_lines_modified_sites",which_nucl,".xls")),sheetName="all_significant_sites",row.names=F,col.names=T)
  
final_id_U <- as.data.frame(bind_rows(final_id)) %>% 
  subset(junction=="No junction") %>% 
  subset(grepl(which_nucl,ref_kmer)==T)%>%
  subset(canonicity=="C")%>%
  dplyr::rename(genomicStartKmer = genomicPos) %>%
  mutate(genomicEndKmer = (genomicStartKmer+5)) %>%
  relocate(genomicEndKmer, .after = genomicStartKmer)                                            # ALL CANONICAL Us : significant,no junctions, canonical, U present
write.xlsx(final_id_U, rdp(paste0("all_cell_lines_modified_sites",which_nucl,".xls")),sheetName=paste0("all_canonical_",which_nucl,"s_significant_sites"),row.names=F,col.names=T,append=TRUE)
  
final_id_5p <- as.data.frame(bind_rows(final_id)) %>% 
  subset(genomicPos<=100)%>%
  dplyr::rename(genomicStartKmer = genomicPos) %>%
  mutate(genomicEndKmer = (genomicStartKmer+5)) %>%
  relocate(genomicEndKmer, .after = genomicStartKmer)                                             # 5p sites : significant,genomicpos <=100
write.xlsx(final_id_5p, rdp(paste0("all_cell_lines_modified_sites",which_nucl,".xls")),sheetName="5p_significant_sites",row.names=F,col.names=T,append=TRUE)
final_id_5p_Us <- as.data.frame(bind_rows(final_id)) %>% 
  subset(genomicPos<=100) %>%
  subset(canonicity=="C") %>%
  subset(grepl(which_nucl,ref_kmer)==T)%>%
  dplyr::rename(genomicStartKmer = genomicPos) %>%
  mutate(genomicEndKmer = (genomicStartKmer+5)) %>%
  relocate(genomicEndKmer, .after = genomicStartKmer)                                             # 5p Us : significant,U present, canonical, genomicpos <=100
write.xlsx(final_id_5p_Us, rdp(paste0("all_cell_lines_modified_sites",which_nucl,".xls")),sheetName=paste0("5p_canonical_",which_nucl,"s_sign_sites"),row.names=F,col.names=T,append=TRUE)
  
# Table with Fleming sites
final_fleming <- lapply(toplot,function(x){
  x <- x %>%
    subset(fleming_presence==T)%>%
    dplyr::rename(genomicStartKmer = genomicPos) %>%
    mutate(genomicEndKmer = (genomicStartKmer+5)) %>%
    relocate(genomicEndKmer, .after = genomicStartKmer) %>%
    dplyr::select(-SAMPLEID,-start_end_sites,-modif,-others)
  return(x)
})
final_fleming <- as.data.frame(bind_rows(final_fleming))                        # Fleming ALL sites 
write.xlsx(final_fleming, rdp(paste0("all_cell_lines_modified_sites",which_nucl,".xls")),sheetName="Fleming_redundant_ALL_sites",row.names=F,col.names=T,append=TRUE)
  
final_fleming_sign <- lapply(toplot,function(x){
  x <- x %>%
    subset(fleming_presence==T & significant==T)%>%
    dplyr::rename(genomicStartKmer = genomicPos) %>%
    mutate(genomicEndKmer = (genomicStartKmer+5)) %>%
    relocate(genomicEndKmer, .after = genomicStartKmer) %>%
    dplyr::select(-SAMPLEID,-start_end_sites,-modif,-others)                         # Fleming significant sites 
  return(x)
})
final_fleming_sign <- as.data.frame(bind_rows(final_fleming_sign))
write.xlsx(final_fleming_sign, rdp(paste0("all_cell_lines_modified_sites",which_nucl,".xls")),sheetName="Fleming_redundant_sign_sites",row.names=F,col.names=T,append=TRUE)
  
final_fleming_unique <- final_fleming_sign %>%
  dplyr::select(genomicStartKmer,genomicEndKmer,ref_kmer)%>%
  dplyr::group_by(genomicStartKmer,genomicEndKmer,ref_kmer) %>% 
  dplyr::filter(row_number() == 1)
final_fleming_unique<-as.data.frame(final_fleming_unique)                       # Fleming significant, non redundant sites
write.xlsx(final_fleming_unique, rdp(paste0("all_cell_lines_modified_sites",which_nucl,".xls")),sheetName="Fleming_non_red_CandNC_sign_sites",row.names=F,col.names=T,append=TRUE)
  
final_fleming_unique<- final_fleming %>%
  dplyr::select(genomicStartKmer,genomicEndKmer,ref_kmer)%>%
  dplyr::group_by(genomicStartKmer,genomicEndKmer,ref_kmer) %>% 
  dplyr::filter(row_number() == 1)
final_fleming_unique<-as.data.frame(final_fleming_unique)                       # Fleming non redundant, ALL sites
write.xlsx(final_fleming_unique, rdp(paste0("all_cell_lines_modified_sites",which_nucl,".xls")),sheetName="Fleming_non_red_CandNC_all_sites",row.names=F,col.names=T,append=TRUE)
  
