library(tidyverse)
library(tidyr)
library(stringr)
library(dplyr)
library(xlsx)
library(seqinr)
library(rstudioapi)
library(readr)

source(paste(dirname(getSourceEditorContext()$path),"variables.R",sep="/"))     # source script with variables
source(paste(dirname(getSourceEditorContext()$path),"functions.R",sep="/"))     # source script with functions

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


fragments_bed <- read_tsv(fragments_bedfile, col_types = "cnn") %>%             # load genomic fragments bedfile
  dplyr::mutate(start = (start - 1), end = (end - 1))

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
  mutate(fragment_ID = get_fragment_partial_overlap(genomicPos, fragments_bed)) %>%  # check if any nucleotide of the k-mer overlaps a genomic fragment
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

# overlap of gRNAs k-mers with sgRNAs shared sites identified previously
sgrnas_files <- list.files(SGRNASDIR, pattern="3_shared*", full.names = TRUE) 
sgrna_sites <- lapply(sgrnas_files,function(x){
  x<- read_tsv(x,col_names = c("chr","start","end","name","score","strand"))
  return(x)
}) 
sgrna_sites <- bind_rows(sgrna_sites)

grna_sites <- total %>%
  subset(GMM_logit_pvalue<=pval_thresh_gRNAs & junctions=="No_junction" & abs(Logit_LOR)>=LOR_thresh) %>%
  mutate(sgintersection=kmer_overlaps_sgRNA_site(genomicPos,sgrna_sites))



##############################  PLOTS ##############################
pdf(rdp("WT_vs_IVT_sharkfin_no_junctions.pdf"),height=15,width=20)
  total %>%
    subset(junctions=="No_junction") %>%
    {
      ggplot(., aes(x=abs(Logit_LOR), y=-log10(GMM_logit_pvalue),color=fragment_ID)) +
        geom_point(size=2) +
        {if (nrow(subset(total,(GMM_logit_pvalue<=pval_thresh_gRNAs & abs(Logit_LOR)>=LOR_thresh)))>0) ggrepel::geom_label_repel(data=dplyr::filter(., (GMM_logit_pvalue<=pval_thresh_gRNAs & abs(Logit_LOR)>=LOR_thresh & grepl("U",ref_kmer)==T)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=3,max.overlaps = Inf,force = 8)}+
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
write.table(total,file=rdp("WT_vs_IVT_nanocompore.tsv"),sep = "\t", quote = F, col.names = T, row.names = F)
write.table(grna_sites,file=rdp("WT_vs_IVT_sign.tsv"),sep = "\t", quote = F, col.names = T, row.names = F)


############## BED FILES FOR UCSC GENOME BROWSER TRACK ##############
grna_bed <- grna_sites %>% 
  dplyr::select(ref_id,genomicPos,ref_kmer)%>%
  mutate(ref_id="NC_045512.2",start=(genomicPos),end=(genomicPos+5),id=paste0("ref_kmer=",ref_kmer),score=0,strand="+")%>%
  dplyr::select(-ref_kmer,-genomicPos)
  
write.table(grna_bed,file=rdp("WT_vs_IVT_significant_gRNAs_nojunctions.bed"),sep = "\t", quote = F, col.names = F, row.names = F)

full_genome_bed <- grna_bed %>%
  mutate(ref_id="NC_045512.2",start=1,end=29903,id="genomicRNA",score=0,strand="+")%>%
  head(1)
write.table(full_genome_bed ,file=rdp("viral_gRNA.bed"),sep = "\t", quote = F, col.names = F, row.names = F)
