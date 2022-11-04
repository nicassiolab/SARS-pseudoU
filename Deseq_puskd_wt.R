library(tidyverse)
library(dplyr)
library(ggpubr)
library(tidyr)
library(R.utils)
library(reshape2)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(purrr)
library(pheatmap)
library(rlist)


############################### FUNCTIONS ######################################

# Function the returns full path from basedir
bdp <- function(relpath){
  return(paste0(ROOTDIR,"/",relpath))
}

bdp2 <- function(relpath){
  return(paste0(ROOTDIR2,"/",relpath))
}

rdp <- function(relpath){
  return(paste0(RESULTSDIR,"/",relpath))
}

# Function the returns gene name from the gtf
ext_gene_name <- function(x){
  x<-str_extract(x, '(?<=gene_name\\s)\\w+')
  return(x)
}
ext_gene_id <- function(x){
  x<-str_extract(x, '(?<=gene_id\\s)\\w+')
  return(x)
}
ext_tx_id <- function(x){
  x<-str_extract(x, '(?<=transcript_id\\s)\\w+')
  return(x)
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

############################  PATHS ############################################

## Datasets

# SARS-CoV-2 WT C34 (DRS_CaCo2_3)
# SARS-CoV-2 WT C37 (DRS_CaCo2_4)
# SARS-CoV-2 PUS7KD C34 + doxy (DRS_CaCo2_3)
# SARS-CoV-2 PUS7KD C37 (DRS_CaCo2_4)

KEEP_WT_DRS3=T
KEEP_WT_DRS4=T
KEEP_PUS7KD_DRS3=T
KEEP_PUS7KD_DRS4=T
viral= F
human= T
USE_TPMs=T

ROOTDIR="/Volumes/scratch/FN/TL/cugolini/cov/analysis"
ROOTDIR2="/Volumes/scratch/TSSM/cugolini/cov"
if(viral==T){
  RESULTSDIR=bdp2("analysis/quantification/CaCo2/DESeq/DRS3_DRS4/viral")
  COUNTS_DRS3=bdp2("analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOUNT/counts")
  COUNTS_DRS4=bdp2("analysis/PUS7_KD/map_to_recap_assembly/NANOCOUNT/counts")}
if(human==T){RESULTSDIR=bdp2("analysis/quantification/CaCo2/DESeq/DRS3_DRS4/human")
  COUNTS_DRS3=bdp2("analysis/PUS7_KD_C37/map_to_human_transcriptome/NANOCOUNT/counts")
  COUNTS_DRS4=bdp2("analysis/PUS7_KD/map_to_human_transcriptome/NANOCOUNT/counts")}
HUMAN_GTF="/Volumes/scratch/FN/camilla/nanopore/data/Homo_sapiens.GRCh38.98.gtf"
VIRAL_ASSEMBLY_BED=bdp2("scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus.bed")
ANNOTATED_VIRAL_ASSEMBLY_BED=bdp2("analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/orf_annotate/orf_annotate.bed")

dir.create(RESULTSDIR,recursive = T)

#################################  CODE  #######################################

# Assembly processing
tx <- read_tsv(ANNOTATED_VIRAL_ASSEMBLY_BED, col_names=c("chr", "start", "end", "name", "score", "strand", "cdsStart", "cdsEnd", ".", "ex", "exLen", "exSt"), col_types="cnncncnnnncc")

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

tx_lengths <-read.table(VIRAL_ASSEMBLY_BED, sep = '\t',header = FALSE) %>%
  separate(V11, into=c("ex1","ex2","ex3") ,sep=",", remove=F) 
tx_lengths[is.na(tx_lengths)] <- 0
tx_lengths <- as.data.frame(tx_lengths) %>% 
  mutate(sumrow= as.numeric(ex1)  + as.numeric(ex2)+as.numeric(ex3)) %>%
  dplyr::rename(length=sumrow) %>%
  dplyr::rename(id=V4) %>%
  select(id,length)

assembly <-read.table(VIRAL_ASSEMBLY_BED, col.names = c("chrom","start","end","id","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"), sep="\t")
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



# Human transcriptome annotation processing
gtf <- read.table(HUMAN_GTF, header = FALSE, sep = '\t')
gtf$gene_name <- sapply(gtf$V9,ext_gene_name)
gtf$gene_id <- sapply(gtf$V9,ext_gene_id)
gtf$transcript_id <- sapply(gtf$V9,ext_tx_id)
gtf_tx <- subset(gtf,gtf$V3 == "transcript") %>%
  select(gene_id,transcript_id,gene_name) 



# Load the data
if(USE_TPMs==T){count_to_use="tpm"}else{count_to_use="est_count"}

sampleDRS4_list <- list.files(COUNTS_DRS4, full.names = TRUE)
sampleDRS4_counts <- lapply(sampleDRS4_list, function(x){
  sample_name <- gsub("\\counts.tsv","DRS4",x)
  sample_name <- gsub(".*/","",sample_name)
  x<-read.csv(x,header = TRUE, sep = "\t") %>%
    select(transcript_name,all_of(count_to_use)) %>%
    dplyr::rename(!!sample_name:=all_of(count_to_use))
})

sampleDRS3_list <- list.files(COUNTS_DRS3, full.names = TRUE)
sampleDRS3_counts <- lapply(sampleDRS3_list, function(x){
  sample_name <- gsub("\\counts.tsv","DRS3",x)
  sample_name <- gsub(".*/","",sample_name)
  x<-read.csv(x,header = TRUE, sep = "\t") %>%
    select(transcript_name,all_of(count_to_use)) %>%
    dplyr::rename(!!sample_name:=all_of(count_to_use))
})



####################### CODE  #############################
sample_counts <- list()
a <- NULL
b <- NULL
c <- NULL
WT_samples_in_filename <- c()
KD_samples_in_filename <- c()

if(KEEP_WT_DRS3==T){
  sample_counts <- list.append(sample_counts,sampleDRS3_counts[[2]])
  a <- append(a,"WT_DRS3")
  b <- append(b,'wt')
  c <- append(c,"CaCo2")
  WT_samples_in_filename <- c("DRS3",WT_samples_in_filename)
}
if(KEEP_WT_DRS4==T){
  sample_counts <- list.append(sample_counts,sampleDRS4_counts[[2]])
  a <- append(a,"WT_DRS4")
  b <- append(b,'wt')
  c <- append(c,"CaCo2")
  WT_samples_in_filename <- c("DRS4",WT_samples_in_filename)
}
if(KEEP_PUS7KD_DRS3==T){
  sample_counts <- list.append(sample_counts,sampleDRS3_counts[[1]])
  a <- append(a,"PUS7KD_DRS3")
  b <- append(b,'kd')
  c <- append(c,"CaCo2")
  KD_samples_in_filename <- c("DRS3",KD_samples_in_filename)
}
if(KEEP_PUS7KD_DRS4==T){
  sample_counts <- list.append(sample_counts,sampleDRS4_counts[[1]])
  a <- append(a,"PUS7KD_DRS4")
  b <- append(b,'kd')
  c <- append(c,"CaCo2")
  KD_samples_in_filename <- c("DRS4",KD_samples_in_filename)
}


# Create metadata
metaData <- data.frame(a,b,c)         # THE ORDER OF THE METADATA MUST BE THE SAME OF THE TABLE
colnames(metaData)<- c("id","dex","celltype")
# Create counts dataset
countData <- sample_counts %>%
  purrr::reduce(full_join, by="transcript_name")
countData[is.na(countData)] <- 0
countData[2:ncol(countData)] <- lapply(countData[2:ncol(countData)], as.integer)


if(USE_TPMs==F){
  # Transform data
  # create object
  dds <- DESeqDataSetFromMatrix(countData=countData, 
                                colData=metaData, 
                                design=~dex, tidy = TRUE)
  
  # Run Deseq2
  # run DESeq
  dds <- DESeq(dds)
  # display DESeq results 
  res <- results(dds)
  
  # Save results to file
  # sort results by pvalue
  resOrdered <- res[order(res$pvalue),]
  
  # write file
  write.table(
    as.data.frame(resOrdered),
    file = rdp(
      paste0(
        "WT_",
        paste0(WT_samples_in_filename, collapse = "_"),
        "_PUS7KD_",
        paste0(KD_samples_in_filename, collapse = "_"),
        "_results.csv"
      )
    ),
    quote = F,
    sep = '\t'
  )
  
  
  # Data transformation
  # data transformation
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  rld <- rlog(dds, blind=FALSE)
  # calculate sample distances
  sampleDists <- dist(t(assay(vsd)))
  # set data for Sample Distance Matrix
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  ############################  PLOTS ############################################
  
  pdf(file=rdp(
    paste0(
      "WT_",
      paste0(WT_samples_in_filename, collapse = "_"),
      "_PUS7KD_",
      paste0(KD_samples_in_filename, collapse = "_"),
      "_plots.pdf"
    )
  ))
  
  #  Volcano Plot
  par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
  topT <- as.data.frame(res)
  #Adjusted P values (FDR Q values)
  with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
  with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
  #with(subset(topT, padj<1 ), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<1 ), cex=0.8, pos=3))
  #Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
  abline(v=0, col="black", lty=3, lwd=1.0)
  abline(v=-2, col="black", lty=4, lwd=2.0)
  abline(v=2, col="black", lty=4, lwd=2.0)
  #abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
  
  plotPCA(vsd, intgroup=c("dex"))   #PCA
  plotDispEsts(dds)                 #dispersion plot
  ##questo plot è usato per vedeere come si distribuisce la varianza rispetto all'espressione. 
  ##la relazione sottesa da DeSeq è che all'aumentare dell'espressione diminuisce la varianza 
  ##(in questo caso la dispersion=varianza/media) perciò la curva fittata dev'essere un'iperbole. 
  ##I punti in nero sono i punti dei nostri dati , mentre la curva rossa è quella che li fitta e 
  ##i punti blu sono gli outliers per i quali il metodo statistico di Deseq decide di assegnare 
  ##il valore degli estimated invece che quello del fit. Infatti normalmente Deseq approssima 
  ##la varianza biologica e il rumore sotteso tramite il fit assegnando ai punti la coordinata del 
  ##fit e non quella degli estimated. Perciò questa distribuzione deve avere sempre la forma di 
  ##un'iperbole o significa che qualcosa è errato
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)              #Sample Distance Matrix
  
  dev.off()
}





# As samples are extremely different one from the other, we try with TPMs 
# instead of est_count, thus avoiding DESeq internal normalisation

countData <- countData %>% 
  mutate(DRS3fc=PUS7_KD_DRS3/WT_DRS3)%>% 
  mutate(DRS4fc=PUS7_KD_DRS4/WT_DRS4)%>%
  dplyr::rename(id=transcript_name)

if(viral==T){
  
  canonicity <- canonicity%>%
    dplyr::rename(id=ref_id)
  
  countData <- left_join(countData,names, by="id")%>%
    separate(id,into = c("id","coord"),sep = "::")%>%
    left_join(canonicity, by="id") %>%
    left_join(tx_lengths, by="id")
  
  canonical_counts <- countData%>%
    subset(canonicity=="C") %>%
    subset(name!="ORF10_SARS2") %>%
    subset(name!="ORF9D_SARS2")
  
  mean_lengths <- canonical_counts%>%
    dplyr::group_by(name)%>%
    dplyr::summarise(length=mean(length)) 
  
  orf_counts <- canonical_counts %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(
      PUS7_KD_DRS3 = sum(PUS7_KD_DRS3),
      PUS7_KD_DRS4 = sum(PUS7_KD_DRS4),
      WT_DRS3 = sum(WT_DRS3),
      WT_DRS4 = sum(WT_DRS4)
    ) %>%
    mutate(DRS3fc = PUS7_KD_DRS3 / WT_DRS3) %>%
    mutate(DRS4fc = PUS7_KD_DRS4 / WT_DRS4)
  
  orf_counts <- orf_counts %>% 
    left_join(mean_lengths,by="name")
  
  #corr between tpms of the two samples in orf counts
  p1 <- ggplot(orf_counts, aes(x=WT_DRS4, y=PUS7_KD_DRS4)) +
    geom_point() +
    geom_smooth(method = "lm")+
    ggrepel::geom_label_repel(data=orf_counts, aes(label=name), colour="black", size=4)+
    theme_bw(18)+
    stat_cor(method = "pearson") +
    scale_y_continuous(label=fancy_scientific)+
    scale_x_continuous(label=fancy_scientific)+
    ggtitle("Scatterplot of counts in DRS4 samples",subtitle = "Points are orfs")
  
  p2 <- ggplot(orf_counts, aes(x=WT_DRS3, y=PUS7_KD_DRS3)) +
    geom_point() +
    geom_smooth(method = "lm")+
    ggrepel::geom_label_repel(data=orf_counts, aes(label=name), colour="black", size=4)+
    theme_bw(18)+
    stat_cor(method = "pearson") +
    scale_y_continuous(label=fancy_scientific)+
    scale_x_continuous(label=fancy_scientific)+
    ggtitle("Scatterplot of counts in DRS3 samples",subtitle = "Points are orfs")
  
  
  #corr between tpms of the two samples in tx counts
  p3<- ggplot(canonical_counts, aes(x=WT_DRS4, y=PUS7_KD_DRS4)) +
    geom_point() +
    geom_smooth(method = "lm")+
    ggrepel::geom_label_repel(data=canonical_counts, aes(label=name), colour="black", size=4)+
    theme_bw(18)+
    stat_cor(method = "pearson") +
    scale_y_continuous(label=fancy_scientific)+
    scale_x_continuous(label=fancy_scientific)+
    ggtitle("Scatterplot of counts in DRS4 samples",subtitle = "Points are tx counts")
  
  p4 <- ggplot(canonical_counts, aes(x=WT_DRS3, y=PUS7_KD_DRS3)) +
    geom_point() +
    geom_smooth(method = "lm")+
    ggrepel::geom_label_repel(data=canonical_counts, aes(label=name), colour="black", size=4)+
    theme_bw(18)+
    stat_cor(method = "pearson") +
    scale_y_continuous(label=fancy_scientific)+
    scale_x_continuous(label=fancy_scientific)+
    ggtitle("Scatterplot of counts in DRS3 samples",subtitle = "Points are tx counts")
  
  
  #corr between tpms and length of the two samples in tx counts
  p5<-ggplot(orf_counts, aes(x=WT_DRS3, y=1/length)) +
    geom_point() +
    geom_smooth(method = "lm")+
    ggrepel::geom_label_repel(data=orf_counts, aes(label=name), colour="black", size=4)+
    theme_bw(18)+
    stat_cor(method = "pearson") +
    scale_y_continuous(label=fancy_scientific)+
    scale_x_continuous(label=fancy_scientific)+
    ggtitle("Scatterplot of counts in WT DRS3 sample vs 1/length",subtitle = "Points are orf counts")
  
  p6<-ggplot(orf_counts, aes(x=WT_DRS4, y=1/length)) +
    geom_point() +
    geom_smooth(method = "lm")+
    ggrepel::geom_label_repel(data=orf_counts, aes(label=name), colour="black", size=4)+
    theme_bw(18)+
    stat_cor(method = "pearson") +
    scale_y_continuous(label=fancy_scientific)+
    scale_x_continuous(label=fancy_scientific)+
    ggtitle("Scatterplot of counts in WT DRS4 sample vs 1/length",subtitle = "Points are orf counts")
  
  p7<-ggplot(orf_counts, aes(x=PUS7_KD_DRS3, y=1/length)) +
    geom_point() +
    geom_smooth(method = "lm")+
    ggrepel::geom_label_repel(data=orf_counts, aes(label=name), colour="black", size=4)+
    theme_bw(18)+
    stat_cor(method = "pearson") +
    scale_y_continuous(label=fancy_scientific)+
    scale_x_continuous(label=fancy_scientific)+
    ggtitle("Scatterplot of counts in PUS7KD DRS3 sample vs 1/length",subtitle = "Points are orf counts")
  
  p8<-ggplot(orf_counts, aes(x=PUS7_KD_DRS4, y=1/length)) +
    geom_point() +
    geom_smooth(method = "lm")+
    ggrepel::geom_label_repel(data=orf_counts, aes(label=name), colour="black", size=4)+
    theme_bw(18)+
    stat_cor(method = "pearson") +
    scale_y_continuous(label=fancy_scientific)+
    scale_x_continuous(label=fancy_scientific)+
    ggtitle("Scatterplot of counts in PUS7KD DRS4 sample vs 1/length",subtitle = "Points are orf counts")
  
  #corr between FC of the two samples in orf and tx counts
  
  p9 <- ggplot(canonical_counts, aes(x=DRS3fc, y=DRS4fc)) +
    geom_point() +
    geom_smooth(method = "lm")+
    ggrepel::geom_label_repel(data=canonical_counts, aes(label=name), colour="black", size=4)+
    theme_bw(18)+
    scale_y_continuous(label=fancy_scientific)+
    scale_x_continuous(label=fancy_scientific)+
    stat_cor(method = "pearson", label.x = 0.5)+
    ggtitle("Scatterplot of FC(PUS7KD/WT) in DRS3 and DRS4 samples",subtitle = "Points are canonical transcript models")
  
  p10 <- ggplot(orf_counts, aes(x=DRS3fc, y=DRS4fc)) +
    geom_point() +
    geom_smooth(method = "lm")+
    ggrepel::geom_label_repel(data=orf_counts, aes(label=name), colour="black", size=4)+
    theme_bw(18)+
    scale_y_continuous(label=fancy_scientific)+
    scale_x_continuous(label=fancy_scientific)+
    ggtitle("Scatterplot of FC(PUS7KD/WT) in DRS3 and DRS4 samples",subtitle = "Points are orfs") +
    stat_cor(method = "pearson", label.x = 0.5)
  
  
  
  pdf(file=rdp(
    paste0(
      "TPMS_WT_",
      paste0(WT_samples_in_filename, collapse = "_"),
      "_PUS7KD_",
      paste0(KD_samples_in_filename, collapse = "_"),
      "_plots.pdf"
    )
  ),width = 20,height = 10)
  
  ggarrange(p1,p2)
  ggarrange(p3,p4)
  ggarrange(p5,p6,p7,p8,ncol = 2,nrow = 2)
  ggarrange(p9,p10)
  dev.off()

} 

if(human==T){
  countData <- countData %>%
    separate(id,into = c("transcript_id","strand"),sep = "\\(")%>%
    left_join(gtf_tx, by="transcript_id")
  
  pdf(file=rdp(
    paste0(
      "TPMS_WT_",
      paste0(WT_samples_in_filename, collapse = "_"),
      "_PUS7KD_",
      paste0(KD_samples_in_filename, collapse = "_"),
      "_plots.pdf"
    )
  ),width = 20,height = 10)
  
  
  #corr between FC of the two samples in orf and tx counts
  countData %>%
    subset(DRS3fc!="NaN")%>%
    subset(DRS4fc!="NaN")%>%
    ggplot(., aes(x=DRS3fc, y=DRS4fc)) +
    geom_point() +
    theme_bw()+
    scale_y_continuous(label=fancy_scientific)+
    scale_x_continuous(label=fancy_scientific)+
    ggtitle("Scatterplot of FC(PUS7KD/WT) in DRS3 and DRS4 samples")
  
  dev.off()
  
  # write file
  write_excel_csv2(
    countData,
    file = rdp(
      paste0(
        "TPMS_WT_",
        paste0(WT_samples_in_filename, collapse = "_"),
        "_PUS7KD_",
        paste0(KD_samples_in_filename, collapse = "_"),
        "_results.csv"
      )
    ),
    ,col_names = T
  )

}



