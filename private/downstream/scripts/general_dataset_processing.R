
########################### GENERAL DATA #######################################

# IVT bedfile to use for IVT junctions
IVT <- read_tsv(IVT_bedfile, col_names = c("chrom","start","end","id","score","strand"),col_types="cnncnc") %>% 
  dplyr::mutate(start=(start-1),end=(end-1))

# NRCeq assembly bedfile 
assembly <- read.table(assembly_bedfile, col.names = c("chrom","start","end","id","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"), sep="\t") %>%
  separate(id,into = c("id","others"),sep="::") %>%
  dplyr::select(-others)

# datasets to calculate the length of NRCeq reference sgRNAs
tx_lengths <- assembly %>% 
  mutate(blockSizes = substring(blockSizes,1, nchar(blockSizes)-1))%>%
  separate(blockSizes, into=c("ex1","ex2","ex3") ,sep=",", remove=F,)%>%
  as.data.frame() 
tx_lengths[is.na(tx_lengths)] <- 0
tx_lengths <- tx_lengths %>%
  mutate(sumrow= as.numeric(ex1) + as.numeric(ex2)+as.numeric(ex3))

# ORF annotation bedfile
tx <- read_tsv(ORF_annotation_bedfile, col_names=c("chr", "start", "end", "ref_id", "score", "strand", "cdsStart", "cdsEnd", ".", "ex", "exLen", "exSt"), col_types="cnncncnnnncc")

# genomic fragments bedfile
fragments <- read_tsv(fragments_bedfile,col_types = "cnn")
fragments_bed <- fragments %>% 
  dplyr::mutate(start=(start-1),end=(end-1))

# reference viral genome fasta file
viral_genome <- read.fasta(viral_genome_fasta) %>%
  getSequence(as.string=T) %>% unlist()


# sitelist and motifs according to the modified nucleotide
if(which_nucl=="U"){
  sitelist <- read_tsv(pseudoU_motifs,col_types="ccc") %>%
    dplyr::rename(motif=IUPAC) %>%
    subset(motif=="UNUAR")                                                      # use only UNUAR motifs
}


############################# PROCESS GENERAL DATASETS #########################

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
names <- names %>% 
  dplyr::select(ref_id,ORF) %>% 
  separate(ref_id,into=c("id","others"),sep = "::",remove=T)

assembly <- left_join(assembly,canonicity,by="id")
assembly <- left_join(assembly,names,by="id")

assembly_multiple_isoform <- assembly %>% 
  subset(canonicity=="C") %>% 
  group_by(ORF) %>% 
  filter(n()>1)

