#!/bin/python
# usage: singularity exec -B /hpcnfs/scratch/ /hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_v1.0.4.sif python3 pseudoU_signal_plots.py
import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
#from Bio import SeqIO
from nanocompore.SampCompDB import SampCompDB

transcriptome="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD_WT/'
resdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD_WT/plots/significant_UNUAR/'


#genmodpos=29814
#modpos=1569
#kmer='AUUAA'
#tx='31f63690-13b9-4e8e-bca5-6a5019782682_1194_NC_045512v2:28245_29853'
#tx1='31f63690-13b9-4e8e-bca5-6a5019782682|1194::NC_045512v2:28245-29853'
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos)+'_'+kmer+'.pdf')
#s=SampCompDB(rootdir+"/outSampComp.db", fasta_fn=transcriptome)
#p1 = s.plot_signal(ref_id=tx1, start=(modpos-15), end=(modpos+15), plot_style='seaborn-whitegrid', figsize=[18,9])
#pdf.savefig(p1[0])
#pdf.close()


#genmodpos=29657
#modpos=1412
#kmer='ACUUU'
#tx='31f63690-13b9-4e8e-bca5-6a5019782682_1194_NC_045512v2:28245_29853'
#tx1='31f63690-13b9-4e8e-bca5-6a5019782682|1194::NC_045512v2:28245-29853'
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos)+'_'+kmer+'.pdf')
#s=SampCompDB(rootdir+"/outSampComp.db", fasta_fn=transcriptome)
#p1 = s.plot_signal(ref_id=tx1, start=(modpos-15), end=(modpos+15), plot_style='seaborn-whitegrid', figsize=[18,9])
#pdf.savefig(p1[0])
#pdf.close()

genmodpos=49
modpos=35
kmer='CUUGU'
tx='e94116d8-72eb-4252-819c-ad06fd9870d4_4145_NC_045512v2_14-29871'
tx1='e94116d8-72eb-4252-819c-ad06fd9870d4|4145::NC_045512v2:14-29871'
pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos)+'_'+kmer+'.pdf')
s=SampCompDB(rootdir+"/outSampComp.db", fasta_fn=transcriptome)
p1 = s.plot_signal(ref_id=tx1, start=(modpos-15), end=(modpos+15), plot_style='seaborn-whitegrid', figsize=[18,9])
pdf.savefig(p1[0])
pdf.close()





