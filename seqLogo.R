##  script to plot a sequence logo for modification sites in the Burrows et al paper
##  this in order to identify a motif in the sequences
library(stringr)
library(grid)
library(seqLogo)
library(seqinr)
library(Biostrings)

ROOTDIR="/Volumes/scratch/TSSM/cugolini/cov"

# Function the returns full path from basedir
bdp <- function(relpath){
  return(paste0(ROOTDIR,"/",relpath))
}



genome_fa <-
  read.fasta(
    '/Volumes/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa',
    as.string = T,
    set.attributes = FALSE
  )[[1]]

modif_sites <- c(22322,23317,27164,28417,28759,28927,29418) #modif sites takemn from the paper
modif_seq <- vector()
for (i in modif_sites){   #extract the sequence 10 nt upstream and downstream of the modifications
  temp <- substr(genome_fa,i-10,i+10)
  modif_seq <- c(modif_seq,temp)
}

string_set <- DNAStringSet(modif_seq)
PSSM <- consensusMatrix(string_set, as.prob = TRUE)[1:4,]   #count nucleotide frequencies per position
pwm <- makePWM(PSSM)
seqLogo(pwm)
