library(tidyverse)
library(tidyr)
library(stringr)
library(dplyr)
library(xlsx)
library(seqinr)
library(rstudioapi)    

source(paste(dirname(getSourceEditorContext()$path),"functions.R",sep="/"))

ROOTDIR="/Volumes/scratch/TSSM/cugolini/cov"
RESULTSDIR=bdp("analysis/sgRNAs_mods_detection/guppy_v601/CaCo2/nanocompore/downstream")
dir.create(RESULTSDIR,showWarnings = F)
NANOCOMPORE_SGRNAS=bdp("analysis/sgRNAs_mods_detection/guppy_v601/CaCo2/nanocompore/alignments_to_assembly/WT_vs_PUS7KD/outnanocompore_results.tsv")

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
alignment_to <- "assembly"                                                      # possibilities are assembly and human
against_IVT <- F


########################### GENERAL DATA #######################################

IVT <- read_tsv(bdp("scripts_new/backupped_data/IVT_junctions.bed"), col_names = c("chrom","start","end","id","score","strand"),col_types="cnncnc") %>% 
  dplyr::mutate(start=(start-1),end=(end-1))
assembly <- read.table(bdp("scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus.bed"), col.names = c("chrom","start","end","id","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"), sep="\t")
fragments <- read_tsv(bdp("scripts_new/backupped_data/RAPID/fragments_genomic_coord_UCSC.txt"),col_types = "cnn")
fragments_bed <- fragments %>% 
  dplyr::mutate(start=(start-1),end=(end-1))

# Dataframe that reports the assembly junction sites
assembly_junction_sites <- assembly %>% 
  dplyr::select(start,end,id,blockSizes,blockStarts) %>% 
  separate(blockSizes, into = c("blSize1","blSize2","blSize3"), sep = ",",convert=T) %>%
  separate(blockStarts, into = c("blStart1","blStart2","blStart3"), sep = ",",convert=T) %>%
  mutate(blStart1 = (blStart1 + start),blStart2 = (blStart2 + start),blStart3 = (blStart3 + start)) %>%
  mutate(blEnd1 = (blStart1 + blSize1),blEnd2 = (blStart2 + blSize2),blEnd3 = (blStart3 + blSize3))

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

sgRNAs <- read_tsv(NANOCOMPORE_SGRNAS) %>%
    select(-genomicPos,-chr, -strand)%>%
    separate(cluster_counts,
             into = c("SAMPLEID", "Y", "Z"),
             sep = "_(?=[0-9])",
             remove = F
    )%>%
    dplyr::select(-Y,-Z)%>%
    subset(SAMPLEID != "NC")%>%
    mutate(ref_kmer = gsub("T", "U", ref_kmer))

total <- sgRNAs %>%
  dplyr::rename(genomicPos=pos) %>%
  mutate(Logit_LOR=as.numeric(Logit_LOR))%>%
  left_join(sitelist, by = "ref_kmer") %>%
  mutate(modif=ifelse(is.na(modif),"Others",modif)) %>%
  separate(ref_id,into = c("id","others"),sep="::") %>%
  left_join(names, by=c("id","others"))

if(alignment_to=="assembly"){
  total <- total %>%
    rowwise() %>%
    mutate(IVT=kmer_overlaps_ORFjunc(genomicPos,id)) %>%
    rowwise() %>%
    mutate(fragment_ID = get_fragment_partial_overlap(genomicPos,(genomicPos+4),fragments_bed))
  total$fragment_ID<-factor(total$fragment_ID,levels=c("No_Fragment","Fragment1","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10"))
}
  
if(against_IVT==T){
  total <- total %>%
    rowwise() %>%
    mutate(IVT = kmer_overlaps_IVT(genomicPos,IVT,IVT_junc_interval_left,IVT_junc_interval_right))
}

if(alignment_to=="assembly"){
  total_split <- split(total,total$id)
}

##############################  PLOTS ##############################
pdf(rdp("WT_vs_PUS7KD_sharkfin.pdf"),height=15,width=20)

lapply(total_split, function(x){
  x %>%
    {
      ggplot(., aes(x=abs(Logit_LOR), y=-log10(GMM_logit_pvalue),color=fragment_ID,shape=IVT)) +
        geom_point(size=2) +
        #{if (nrow(subset(x,(GMM_logit_pvalue<=pval_thresh & abs(Logit_LOR)>=LOR_thresh)))>0) ggrepel::geom_label_repel(data=filter(., (GMM_logit_pvalue<=pval_thresh & abs(Logit_LOR)>=LOR_thresh & grepl("U",ref_kmer)==T)) ,aes(label=paste0(ref_kmer, " (",genomicPos,")")), colour="black", size=2)}+
        scale_color_manual(breaks = c("No_Fragment","Fragment1","Fragment1_2","Fragment2_3","Fragment3_4","Fragment4_5","Fragment5","Fragment6","Fragment6_7","Fragment7_8","Fragment8_9","Fragment9_10","Fragment10"),values=c("black","blue", "green","grey","gold","coral","aquamarine","darkgreen","navy","deeppink","magenta","cyan","orange")) +
        scale_shape_manual(values=c(16, 18))+
        theme_bw(22) +
        geom_hline(yintercept=-log10(pval_thresh), linetype="dashed", color = "red")+
        geom_vline(xintercept=LOR_thresh, linetype="dashed", color = "red") +
        ggtitle(label =paste0("Significant k-mers in sgRNA ",unique(x$id)), subtitle = unique(x$ORF))
    }
})
 
dev.off()

##############################  TABLES ##############################
write.table(total,file=rdp("WT_vs_PUS7KD_nanocompore.tsv"),sep = "\t", quote = F, col.names = T, row.names = F)

