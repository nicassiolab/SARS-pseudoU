library(knitr)
library(rmdformats)
library(rmarkdown)
library(tidyverse)
library(ggpubr)
library(GGally)
library(tidyr)
library(R.utils)
library(reshape2)

## Global options
options(max.print="75")
opts_chunk$set(echo=FALSE,
               cache=TRUE,
               prompt=FALSE,
               tidy=FALSE,
               comment=NA,
               message=FALSE,
               warning=FALSE,
               dev=c("png","pdf"))
opts_knit$set(width=75)

ROOTDIR="/Volumes/scratch/TSSM/cugolini/cov"

# Function the returns full path from basedir
bdp <- function(relpath){
  return(paste0(ROOTDIR,"/",relpath))
}


fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}


vlab <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

### GENERAL DATA 

tx <- read_tsv(bdp("analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/orf_annotate/orf_annotate.bed"), col_names=c("chr", "start", "end", "name", "score", "strand", "cdsStart", "cdsEnd", ".", "ex", "exLen", "exSt"), col_types="cnncncnnnncc")

names <- select(tx, orig=name) %>% separate(orig, into=c("id", "protein"), sep="#", remove=F) %>%
  mutate(protein=case_when(is.na(protein)~"Unknown", T~protein)) %>%
  separate(protein, into=c("sp", "uniprot_id", "protein"), sep="\\|", remove=F) %>%
  select(-sp) %>%
  mutate(name=gsub("([^\\(]+).+", "\\1", protein),
         tip=as.numeric(gsub("([^\\(]+)\\(([^%]+)%/([^%]+)%\\)", "\\2", protein)),
         qip=as.numeric(gsub("([^\\(]+)\\(([^%]+)%/([^%]+)%\\)", "\\3", protein))) %>%
  select(-protein) %>%
  select(-orig) %>%
  mutate(name=case_when(is.na(name)~"Unknown", T~name))

names$name[names$id == "efad7b96-ac2e-4ce1-9b83-863ffdb18eac|116::NC_045512v2:11-29873"] <- "ORF10_SARS2"
names$name[names$id == "de81ef19-655d-4ced-a9cc-cb8384001058|107::NC_045512v2:11-29874"] <- "ORF9D_SARS2"

tx_lengths <-read.table(bdp("analysis/recappable_assembly/two_datasets/assemblies/pinfish/aln_consensus.bed"), sep = '\t',header = FALSE) %>%
  separate(V11, into=c("ex1","ex2","ex3") ,sep=",", remove=F) 
tx_lengths[is.na(tx_lengths)] <- 0
tx_lengths <- as.data.frame(tx_lengths) %>% mutate(sumrow= as.numeric(ex1)  + as.numeric(ex2)+as.numeric(ex3))

#### blacklist IVT junctions

sitelist <- read_tsv(bdp("scripts_new/backupped_data/sites_blacklist.txt")) 
colnames(sitelist) <- c("ref_kmer","Modification type","modif")

IVT <-read_tsv(bdp("scripts_new/backupped_data/IVT_junctions.bed"), col_names = c("chrom","start","end","id","score","strand"))
junction_sites <- IVT %>% select(start,end) %>% unlist(use.names=FALSE) %>% sort()
blacklist_IVT <- vector()
for (i in junction_sites){
  temp <- seq(from = i-25, to = i+25, by = 1)
  blacklist_IVT <- c(blacklist_IVT,temp)
}
blacklist_IVT <- blacklist_IVT[which(blacklist_IVT > 0)]

#### blacklist exon junctions

assembly <-read.table(bdp("scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus.bed"), col.names = c("chrom","start","end","id","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"), sep="\t")
assembly_junction_sites <- assembly %>% select(start,end,id,blockSizes,blockStarts) %>% 
  separate(blockSizes, into = c("blSize1","blSize2","blSize3"), sep = ",") %>%
  separate(blockStarts, into = c("blStart1","blStart2","blStart3"), sep = ",") %>%
  mutate(blStart1 = (as.numeric(blStart1) + as.numeric(start))) %>%
  mutate(blStart2 = (as.numeric(blStart2) + as.numeric(start))) %>%
  mutate(blStart3 = (as.numeric(blStart3) + as.numeric(start))) %>%
  mutate(blEnd1 = (as.numeric(blStart1) + as.numeric(blSize1))) %>%
  mutate(blEnd2 = (as.numeric(blStart2) + as.numeric(blSize2))) %>%
  mutate(blEnd3 = (as.numeric(blStart3) + as.numeric(blSize3))) 

#canonicity
canonicity <- select(assembly,start,end,id) %>%
  mutate(canonicity=ifelse(start>100,"NC", "C")) %>%
  mutate(canonicity=ifelse(end<29000,"NC", canonicity)) %>%
  select(-start,-end)
colnames(canonicity) <- c("ref_id","canonicity") 


#### functions

pvalues_fisher_method = function(pvalues){
  keep = (pvalues >= 0) & (pvalues <= 1)
  pvalues[pvalues == 0] = 1e-285
  lnp = log(pvalues)
  chisq = (-2) * rowSums(lnp)
  df = 2 * length(lnp)
  fisher_pval = stats::pchisq(chisq, df, lower.tail=FALSE)
  return(fisher_pval)
}

junc_blacklist = function(tx) {
  corr_row <- subset(assembly_junction_sites, assembly_junction_sites$id %in% tx)
  blacklist_vec <- as.numeric(as.vector(corr_row[1,7:12]))
  blacklist_vec <- blacklist_vec[!is.na(blacklist_vec)]
  black_int <- vector()
  for (i in blacklist_vec){
    temp <- seq(from = i-15, to = i+15, by = 1)
    black_int <- c(black_int,temp)
  }
  black_int <- black_int[which(black_int > 0)] %>% sort()
  return(black_int)
}           


tx_list <- list.files(path = "/Volumes/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_genome/NANOCOMPORE/sampcomp/",pattern = "*_results.tsv" , full.names = TRUE,  recursive = T)

PUS7_KD_WT <- lapply(tx_list, function(x){
  x <- read_tsv(x, col_types="ncccccnnncncn")
  x$SAMPLEID <- (separate(data.frame(A = x$cluster_counts), col = "A" , into = c("X", "Y","Z"), sep = "_(?=[0-9])"))$X
  x <- x[!x$SAMPLEID == "NC", ]
  x <- x %>% mutate(ref_kmer=gsub("T","U", ref_kmer))
})

