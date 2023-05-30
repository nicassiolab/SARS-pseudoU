library(rstudioapi)

source(paste(dirname(getSourceEditorContext()$path),"functions.R",sep="/"))     # source script with functions
source(paste(dirname(getSourceEditorContext()$path),"variables.R",sep="/"))     # source script with variables


dir.create(RESULTSDIR,showWarnings = F)                                         # create results directory


########################### GENERAL DATA #######################################

####   Load datasets for the analysis
IVT <-                                                                          # load Kim et al. IVT bedfile
  read_tsv(
    IVT_bedfile,
    col_names = c("chrom", "start", "end", "id", "score", "strand"),
    col_types = "cnncnc"
  ) %>%
  dplyr::mutate(start = (start - 1), end = (end - 1))                           # use 0-based coordinates

# Fragments file
fragments <- read_tsv(fragments_bedfile) %>%
  mutate(start=(start-1),end=(end-1))


SNPs <- read_tsv(file = SNPs_file) %>%                                          # load file containing known viral strains SNPs
  arrange(genomicPos) %>%
  mutate(genomicPos=(genomicPos-1)) %>%
  dplyr::select(genomicPos) %>%
  unique()


# Load modification motifs depending on the nucleotide 
if(which_nucl=="U") {
  sitelist <- read_tsv(pseudoU_motifs, col_types = "ccc") %>%
    dplyr::rename(motif = IUPAC) %>%
    subset(motif == "UNUAR")                                                    # use only UNUAR motifs
} else if (which_nucl == "A") {
  sitelist <- read_tsv(sitelist_m6A_file, col_types = "ccc") %>%
    dplyr::rename(motif = IUPAC) %>%                                          
    subset(modif == "m6A")
}

####################### PROCESSING OF THE DATABASES ############################

# Process gRNAs reads (WT vs IVT)
gRNAs <- read_tsv(nanocompore_gRNAs_file) %>%
    dplyr::select(-genomicPos,-chr, -strand)%>%
    separate(cluster_counts,
             into = c("SAMPLEID", "Y", "Z"),
             sep = "_(?=[0-9])",
             remove = F
    )%>%
    dplyr::select(-Y,-Z)%>%
    subset(SAMPLEID != "NC")%>%
    mutate(ref_kmer = gsub("T", "U", ref_kmer))

total <- gRNAs %>%
  dplyr::rename(genomicPos = pos) %>%
  mutate(Logit_LOR = as.numeric(Logit_LOR)) %>%
  left_join(sitelist, by = "ref_kmer") %>%
  mutate(modif = ifelse(is.na(modif), "Others", modif)) %>%
  rowwise() %>%
  mutate(junctions = kmer_overlaps_IVT(                                         # check if any nucleotide of the k-mer overlaps an IVT junction site
    genomicPos,
    IVT,
    IVT_junc_interval_left,
    IVT_junc_interval_right
  )) %>%
  rowwise() %>%
  mutate(junctions = ifelse(                                                    # check if any nucleotide of the k-mer overlaps an annotated SNP
    is.null(kmer_overlaps_SNPs(genomicPos, SNPs)),
    junctions,
    kmer_overlaps_SNPs(genomicPos, SNPs)
  )) %>%
  rowwise() %>%
  mutate(fragment_ID = get_fragment_partial_overlap(genomicPos,fragments)) %>%  # check if any nucleotide of the k-mer overlaps a genomic fragment
  rowwise() %>%
  mutate(fleming= get_fleming(genomicPos))
   
total$fragment_ID <-
  factor(
    total$fragment_ID,
    levels = c(
      "No_Fragment",
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
      "Fragment10"
    )
  )

gRNAs_sites <- total %>%
  subset(GMM_logit_pvalue<=pval_thresh_gRNAs & junctions=="No_junction" & abs(Logit_LOR)>=LOR_thresh) 
gRNAs_final <- gRNAs_sites %>%
  dplyr::rename(genomicStartKmer = genomicPos, chr = ref_id) %>%
  mutate(genomicEndKmer = (genomicStartKmer+5)) %>%
  relocate(genomicEndKmer, .after = genomicStartKmer) 
gRNAs_final_U <- gRNAs_final %>%
  subset(grepl("U",ref_kmer)==T)



# overlap of all cell lines pulled sgRNAs with gRNAs
sgRNAs_final <- read.xlsx(SGRNAS_all_cell_lines,sheetName="all_canonical_Us_significant_si") %>%
  subset(junction=="No junction" & significant==T & canonicity=="C") 
sgRNAs_final <- unique(sgRNAs_final[c("chr", "genomicStartKmer","genomicEndKmer","ref_kmer","fragment_ID")])
sgRNAs_final_U <- sgRNAs_final %>%
  subset(grepl("U",ref_kmer)==T)


common_sites <- inner_join(gRNAs_final_U,sgRNAs_final_U, by=c("chr","genomicStartKmer","genomicEndKmer","ref_kmer","fragment_ID"))
not_common_sites <- subset(gRNAs_final_U, !(gRNAs_final_U$genomicStartKmer %in% sgRNAs_final_U$genomicStartKmer))

##############################  PLOTS ##############################
pdf(rdp("WT_vs_IVT_sharkfin_no_junctions.pdf"),height=15,width=20)
  total %>%
    subset(junctions=="No_junction") %>%
    {
      ggplot(., aes(x=abs(Logit_LOR), y=-log10(GMM_logit_pvalue),color=fragment_ID)) +
        geom_point(size=2) +
        {if (nrow(subset(total,(GMM_logit_pvalue<=pval_thresh_gRNAs & abs(Logit_LOR)>=LOR_thresh)))>0) ggrepel::geom_label_repel(data=dplyr::filter(., (GMM_logit_pvalue<=pval_thresh_gRNAs & abs(Logit_LOR)>=LOR_thresh & grepl("U",ref_kmer)==T)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=6,max.overlaps = Inf,force = 30)}+
        scale_color_manual(breaks = c("No_Fragment","Fragment1","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
        theme_bw(22) +
        geom_hline(yintercept=-log10(pval_thresh_gRNAs), linetype="dashed", color = "red")+
        geom_vline(xintercept=LOR_thresh, linetype="dashed", color = "red") +
        ggtitle("Significant k-mers (no junctions) in gRNAs")
    }
dev.off()

pdf(rdp("WT_vs_IVT_sharkfin_only_junctions.pdf"),height=15,width=20)
total %>%
  subset(junctions!="No_junction") %>%
  {
    ggplot(., aes(x=abs(Logit_LOR), y=-log10(GMM_logit_pvalue),color=fragment_ID)) +
      geom_point(size=2) +
      {if (nrow(subset(total,(GMM_logit_pvalue<=pval_thresh_gRNAs & abs(Logit_LOR)>=LOR_thresh)))>0) ggrepel::geom_label_repel(data=dplyr::filter(., (GMM_logit_pvalue<=pval_thresh_gRNAs & abs(Logit_LOR)>=LOR_thresh & grepl("U",ref_kmer)==T)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=3)}+
      scale_color_manual(breaks = c("No_Fragment","Fragment1","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
      theme_bw(22) +
      geom_hline(yintercept=-log10(pval_thresh_gRNAs), linetype="dashed", color = "red")+
      geom_vline(xintercept=LOR_thresh, linetype="dashed", color = "red") +
      ggtitle("Significant k-mers (no junctions) in gRNAs")
  }
dev.off()


pdf(rdp("WT_vs_IVT_genomic_plot.pdf"),height=15,width=20)

total %>%
  subset(GMM_logit_pvalue<=pval_thresh_gRNAs & junctions=="No_junction" & abs(Logit_LOR)>=LOR_thresh) %>%
  {
    ggplot(., aes(x=genomicPos, y=-log10(GMM_logit_pvalue)))+
      geom_line() +
      scale_y_continuous(label=fancy_scientific) +
      theme_bw(30) + vlab 
    }
dev.off()

pdf(rdp("WT_vs_IVT_genomic_plot_0_10000.pdf"),height=15,width=20)

total %>%
  subset(GMM_logit_pvalue<=pval_thresh_gRNAs & junctions=="No_junction" & abs(Logit_LOR)>=LOR_thresh) %>%
  {
    ggplot(., aes(x=genomicPos, y=-log10(GMM_logit_pvalue)))+
      geom_line() +
      scale_y_continuous(label=fancy_scientific) +
      theme_bw(30) + vlab +xlim(0,10000)
  }
dev.off()



##############################  TABLES ##############################
total_bed <- total %>%
  dplyr::rename(genomicStartKmer = genomicPos,chr=ref_id) %>%
  mutate(genomicEndKmer = (genomicStartKmer+5)) %>%
  relocate(genomicStartKmer,.after = (chr)) %>%
  relocate(genomicEndKmer,.after = (genomicStartKmer))
write.xlsx(as.data.frame(total_bed), file=rdp("WT_vs_IVT_nanocompore.xls"),
           sheetName="WT_vs_IVT_all_sites",row.names=F,col.names=T)

gRNAs_sites_bed <- gRNAs_sites %>%
  dplyr::rename(genomicStartKmer = genomicPos,chr=ref_id) %>%
  mutate(genomicEndKmer = (genomicStartKmer+5)) %>%
  relocate(genomicStartKmer,.after = (chr)) %>%
  relocate(genomicEndKmer,.after = (genomicStartKmer))
write.xlsx(as.data.frame(gRNAs_sites_bed), file=rdp("WT_vs_IVT_nanocompore.xls"),
           sheetName="WT_vs_IVT_sign_sites",row.names=F,col.names=T,append=T)

gRNAs_final_U <- gRNAs_final_U %>%
  relocate(genomicStartKmer,.after = (chr)) %>%
  relocate(genomicEndKmer,.after = (genomicStartKmer))
write.xlsx(as.data.frame(gRNAs_final_U), file=rdp("WT_vs_IVT_nanocompore.xls"),
           sheetName="WT_vs_IVT_sign_U_sites",row.names=F,col.names=T,append=T)

common_sites <- common_sites %>%
  relocate(genomicStartKmer,.after = (chr)) %>%
  relocate(genomicEndKmer,.after = (genomicStartKmer))
write.xlsx(as.data.frame(common_sites), file=rdp("WT_vs_IVT_nanocompore.xls"),
           sheetName="WT_vs_IVT_sign_U_sites_common_to_sgRNAs",row.names=F,col.names=T,append=T)



############## BED FILES FOR UCSC GENOME BROWSER TRACK ##############
grna_bed <- gRNAs_final_U %>%
  dplyr::select(chr,genomicStartKmer,genomicEndKmer,ref_kmer)%>%
  dplyr::rename(start=genomicStartKmer,end=genomicEndKmer)%>%
  mutate(id=paste0("ref_kmer=",ref_kmer),score=0,strand="+")%>%
  dplyr::select(-ref_kmer)
write.table(grna_bed,file=rdp("WT_vs_IVT_significant_U_gRNAs_nojunctions.bed"),sep = "\t", quote = F, col.names = F, row.names = F)

common_sites_bed <- common_sites %>%
  dplyr::select(chr,genomicStartKmer,genomicEndKmer,ref_kmer)%>%
  dplyr::rename(start=genomicStartKmer,end=genomicEndKmer)%>%
  mutate(id=paste0("ref_kmer=",ref_kmer),score=0,strand="+")%>%
  dplyr::select(-ref_kmer)
write.table(common_sites_bed,file=rdp("WT_vs_IVT_common_significant_U.bed"),sep = "\t", quote = F, col.names = F, row.names = F)


full_genome_bed <- grna_bed %>%
  mutate(ref_id="NC_045512.2",start=1,end=29903,id="genomicRNA",score=0,strand="+")%>%
  dplyr::select(-ref_id)%>%
  head(1)
write.table(full_genome_bed ,file=rdp("viral_gRNA.bed"),sep = "\t", quote = F, col.names = F, row.names = F)


######################  Calculate hypergeometric test ##########################
df <- total %>%
  subset(grepl("U",ref_kmer)==T & junctions=="No_junction" & genomicPos>=SPIKETRSB) 

m <- df %>%
  subset(genomicPos %in% sgRNAs_final_U$genomicStartKmer)%>%
  nrow()

n <- (df %>%
  nrow()) - m

q <- nrow(common_sites)

k <- df %>%
  subset(GMM_logit_pvalue<=pval_thresh_gRNAs & abs(Logit_LOR)>=LOR_thresh) %>%
  nrow()
  

phyper(q=(q-1),m=m,n=n,k=k,lower.tail=F)



