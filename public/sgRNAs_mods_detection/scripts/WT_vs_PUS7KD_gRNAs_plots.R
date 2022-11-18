library(tidyverse)
library(tidyr)
library(stringr)
library(dplyr)
library(xlsx)
library(seqinr)

source(paste(dirname(getSourceEditorContext()$path),"functions.R",sep="/"))

ROOTDIR="/Volumes/scratch/TSSM/cugolini/cov"
RESULTSDIR=bdp("analysis/gRNA_mods/guppy_initial/nanocompore/downstream")
dir.create(RESULTSDIR,showWarnings = F)
NANOCOMPORE_GRNAS=bdp("analysis/gRNA_mods/guppy_initial/nanocompore/all_cell_lines/nanocompore/outnanocompore_results.tsv")

########################### PARAMETERS #########################################

`%!in%` <- Negate(`%in%`)
LOR_thresh <- 0.5
pval_thresh <- 0.01
n_samples <- 1
IVT_junc_interval_left <- 25
IVT_junc_interval_right <- 25
fleming_paper_notation <- c(22322,23317,27164,28417,28759,28927,29418)
fleming_sites_bed <- fleming_paper_notation-3
which_nucl <- "U"                                                               # possibilities are U or A



########################### GENERAL DATA #######################################

IVT <- read_tsv(bdp("scripts_new/backupped_data/IVT_junctions.bed"), col_names = c("chrom","start","end","id","score","strand"),col_types="cnncnc") %>% 
  dplyr::mutate(start=(start-1),end=(end-1))
assembly <- read.table(bdp("scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus.bed"), col.names = c("chrom","start","end","id","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"), sep="\t")
fragments <- read_tsv(bdp("scripts_new/backupped_data/RAPID/fragments_genomic_coord_UCSC.txt"),col_types = "cnn")
fragments_bed <- fragments %>% 
  dplyr::mutate(start=(start-1),end=(end-1))


##  input files that are different depending in which nucleotide the modification is
if(which_nucl=="U"){
  sitelist <- read_tsv(bdp("scripts_new/backupped_data/sites_pseudoU_motifs_blacklist.txt"),col_types="ccc") %>%
    dplyr::rename(motif=IUPAC) %>%
    subset(motif=="UNUAR")                                                      # use only UNUAR motifs
}else if(which_nucl=="A"){
  sitelist <- read_tsv(bdp("scripts_new/backupped_data/sites_blacklist.txt"),col_types="ccc") %>%
    dplyr::rename(motif=IUPAC) %>% subset(modif=="m6A")
}

####################### PROCESSING OF THE DATABASES ############################

gRNAs <- read_tsv(NANOCOMPORE_GRNAS) %>%
    select(-genomicPos,-chr, -strand)%>%
    separate(cluster_counts,
             into = c("SAMPLEID", "Y", "Z"),
             sep = "_(?=[0-9])",
             remove = F
    )%>%
    dplyr::select(-Y,-Z)%>%
    subset(SAMPLEID != "NC")%>%
    mutate(ref_kmer = gsub("T", "U", ref_kmer))

total <- gRNAs %>%
  dplyr::rename(genomicPos=pos) %>%
  mutate(Logit_LOR=as.numeric(Logit_LOR))%>%
  left_join(sitelist, by = "ref_kmer") %>%
  mutate(modif=ifelse(is.na(modif),"Others",modif)) %>%
  rowwise() %>%
  mutate(IVT = kmer_overlaps_IVT(genomicPos,IVT,IVT_junc_interval_left,IVT_junc_interval_right)) %>%
  rowwise() %>%
  mutate(fragment_ID = get_fragment_partial_overlap(genomicPos,(genomicPos+4),fragments_bed))

total$fragment_ID<-factor(total$fragment_ID,levels=c("No_Fragment","Fragment1_5UTR","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10_3UTR"))

##############################  PLOTS ##############################
pdf(rdp("WT_vs_PUS7KD_sharkfin.pdf"),height=15,width=20)
  total %>%
    {
      ggplot(., aes(x=abs(Logit_LOR), y=-log10(GMM_logit_pvalue),color=fragment_ID,shape=IVT)) +
        geom_point(size=2) +
        {if (nrow(subset(total,(GMM_logit_pvalue<=pval_thresh & abs(Logit_LOR)>=LOR_thresh)))>0) ggrepel::geom_label_repel(data=filter(., (GMM_logit_pvalue<=pval_thresh & abs(Logit_LOR)>=LOR_thresh & grepl("U",ref_kmer)==T)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=2)}+
        scale_color_manual(breaks = c("No_Fragment","Fragment1","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
        scale_shape_manual(values=c(17, 16))+
        theme_bw(22) +
        geom_hline(yintercept=-log10(pval_thresh), linetype="dashed", color = "red")+
        geom_vline(xintercept=LOR_thresh, linetype="dashed", color = "red") +
        ggtitle("Significant k-mers (no junctions) in gRNAs")
    }
dev.off()

##############################  TABLES ##############################
write.table(total,file=rdp("WT_vs_PUS7KD_nanocompore.tsv"),sep = "\t", quote = F, col.names = T, row.names = F)

