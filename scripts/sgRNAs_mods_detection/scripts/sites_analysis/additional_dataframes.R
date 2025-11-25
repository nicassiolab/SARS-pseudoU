source(paste(dirname(getSourceEditorContext()$path),"variables.R",sep="/"))

# Dataframe that returns ORFs for every transcript of the assembly
tx <- read_tsv(ORF_annotation_bedfile, col_names=c("chr", "start", "end", "name", 
                                                   "score", "strand", "cdsStart", 
                                                   "cdsEnd", ".", "ex", "exLen", 
                                                   "exSt"), col_types="cnncncnnnncc")
names <- dplyr::select(tx, orig=name) %>% 
  separate(orig, into=c("id", "protein"), sep="#", remove=F) %>%
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
names <- names %>%
  separate(id, into=c("id","others"), sep="::")

# Dataframe to report transcript length
tx_lengths <- tx %>%
  separate(exLen, into=c("ex1","ex2","ex3") ,sep=",", remove=F) %>%
  # dplyr::mutate_all(na_if,"") 
  mutate(ex1=as.numeric(ex1),ex2=as.numeric(ex2),ex3=as.numeric(ex3)) %>%
  replace(is.na(.), 0) %>% 
  mutate(sumrow= ex1 + ex2 + ex3) %>%
  separate(name, into=c("id","others"), sep="::")


# Dataset for the assembly
assembly <- read.table(assembly_bedfile, col.names = c("chrom","start","end",
                                                       "id","score","strand",
                                                       "thickStart","thickEnd",
                                                       "itemRgb","blockCount",
                                                       "blockSizes","blockStarts"), 
                       sep="\t") %>%
  separate(id, into=c("id","others"), sep="::")

# Dataframe that reports canonicity for every transcript of the assembly
canonicity <- dplyr::select(assembly,start,end,id) %>%
  mutate(canonicity=ifelse(start>100,"NC", "C")) %>%
  mutate(canonicity=ifelse(end<29000,"NC", canonicity)) %>%
  mutate(canonicity=ifelse(id=="efad7b96-ac2e-4ce1-9b83-863ffdb18eac|116","C", canonicity)) %>%  #ORF10
  mutate(canonicity=ifelse(id=="de81ef19-655d-4ced-a9cc-cb8384001058|107","C", canonicity)) %>%  #ORF9D
  dplyr::select(-start,-end) 

# Assembly divided per isoform
assembly_named <- assembly %>%
  left_join(names, by=c("id","others")) %>%
  left_join(canonicity, by="id")
assembly_named_canonical <- assembly_named %>%
  subset(canonicity=="C")
n_occur <- data.frame(table(assembly_named_canonical$name))
assembly_multiple_isoform <- assembly_named_canonical %>%
  subset(name %in% as.vector(n_occur[n_occur$Freq > 1,]$Var1))


# Sitelist and motifs for pseudouridines
sitelist <- if(which_nucl == "U"){
  read_tsv(pseudoU_motifs,col_types="ccc") %>%
    dplyr::rename(motif=IUPAC) 
}

# Dataframe that reports the assembly junction sites
assembly_junction_sites <- assembly %>% 
  dplyr::select(start,end,id,blockSizes,blockStarts) %>% 
  separate(blockSizes, into = c("blSize1","blSize2","blSize3"), sep = ",") %>%
  separate(blockStarts, into = c("blStart1","blStart2","blStart3"), sep = ",") %>%
  mutate(blStart1 = (as.numeric(blStart1) + as.numeric(start))) %>%
  mutate(blStart2 = (as.numeric(blStart2) + as.numeric(start))) %>%
  mutate(blStart3 = (as.numeric(blStart3) + as.numeric(start))) %>%
  mutate(blEnd1 = (as.numeric(blStart1) + as.numeric(blSize1))) %>%
  mutate(blEnd2 = (as.numeric(blStart2) + as.numeric(blSize2))) %>%
  mutate(blEnd3 = (as.numeric(blStart3) + as.numeric(blSize3)))


# Dataframe that reports IVT junctions
####   Load datasets for the analysis
IVT <-                                                                          # load Kim et al. IVT bedfile
  read_tsv(
    IVT_bedfile,
    col_names = c("chrom", "start", "end", "id", "score", "strand"),
    col_types = "cnncnc"
  ) %>%
  dplyr::mutate(start = (start - 1), end = (end - 1))                           # use 0-based coordinates

junction_sites <- IVT %>% 
  dplyr::select(start,end) %>% 
  unlist(use.names=FALSE) %>% 
  sort()
blacklist_IVT <- vector()
for (i in junction_sites){
  temp <- seq(from = i-IVT_junc_interval_left, to = i+IVT_junc_interval_right, by = 1)
  blacklist_IVT <- c(blacklist_IVT,temp)
}
blacklist_IVT <- blacklist_IVT[which(blacklist_IVT > 0)]

# Fragments file
fragments <- read_tsv(fragments_bedfile) %>%
  mutate(start=(start-1),end=(end-1))

# Viral genome
viral_genome <- paste(readDNAStringSet(viral_genome_fasta))



