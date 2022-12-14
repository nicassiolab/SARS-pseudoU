#!/bin/python
# usage: plot_results.py
import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
#from Bio import SeqIO
from nanocompore.SampCompDB import SampCompDB

### SIGNALS PLOTS 

transcriptome="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
#fasta_sequences = SeqIO.parse(open(transcriptome),'fasta')
#for fasta in fasta_sequences:
#	print(fasta.id)

## create a pdf for each sample (every sample has n plots = n transcripts)

###############################################################################################################################################################################################################################################################################
											###	PUS7KD VS WT (C37)	###
###############################################################################################################################################################################################################################################################################
### UGUAR
#rootdir = '/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD_WT'
#resdir = '/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD_WT/plots/'
#os.mkdir(resdir)
#kmer='UGUAA'
#s=SampCompDB(rootdir+"/outSampComp.db", fasta_fn=transcriptome)

#genmodpos=51
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')

#with open("/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/corresp_tx_genome_tables/UGUAR_51.txt", "r") as a_file:
#	for line in a_file:
#		split_string = line.split("\t")
#		modpos=int(split_string[0])
#		tx=a_string = split_string[1].rstrip("\n")
#		print(tx)
#		p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,9])
#		pdf.savefig(p1[0])
#pdf.close()

#genmodpos=29079
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')

#with open("/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/corresp_tx_genome_tables/UGUAR_29079.txt", "r") as a_file:
#       for line in a_file:
#               split_string = line.split("\t")
#               modpos=int(split_string[0])
#               tx=a_string = split_string[1].rstrip("\n")
#               print(tx)
#               p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,9])
#               pdf.savefig(p1[0])
#pdf.close()

#genmodpos=29808
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')

#with open("/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/corresp_tx_genome_tables/UGUAR_29808.txt", "r") as a_file:
#       for line in a_file:
#               split_string = line.split("\t")
#               modpos=int(split_string[0])
#               tx=a_string = split_string[1].rstrip("\n")
#               print(tx)
#               p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,9])
#               pdf.savefig(p1[0])
#pdf.close()

#genmodpos=29690
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')

#with open("/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/corresp_tx_genome_tables/UGUAR_29690.txt", "r") as a_file:
#       for line in a_file:
#               split_string = line.split("\t")
#               modpos=int(split_string[0])
#               tx=a_string = split_string[1].rstrip("\n")
#               print(tx)
#               p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,9])
#               pdf.savefig(p1[0])
#pdf.close()

#kmer='UGUAG'
#genmodpos=29649
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')

#with open("/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/corresp_tx_genome_tables/UGUAR_29649.txt", "r") as a_file:
#       for line in a_file:
#               split_string = line.split("\t")
#               modpos=int(split_string[0])
#               tx=a_string = split_string[1].rstrip("\n")
#               print(tx)
#               p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,9])
#               pdf.savefig(p1[0])
#pdf.close()


###	NON-UGUAR MOTIFS

#rootdir = '/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD_WT'
#resdir = '/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD_WT/plots/'
#os.mkdir(resdir)
#kmer='CCAAU'
#s=SampCompDB(rootdir+"/outSampComp.db", fasta_fn=transcriptome)

#genmodpos=24464
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')

#with open("/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/corresp_tx_genome_tables/CCAAU_24464.txt", "r") as a_file:
#       for line in a_file:
#               split_string = line.split("\t")
#               modpos=int(split_string[0])
#               tx=a_string = split_string[1].rstrip("\n")
#               print(tx)
#               p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,9])
#               pdf.savefig(p1[0])
#pdf.close()

#kmer='UAAUA'
#s=SampCompDB(rootdir+"/outSampComp.db", fasta_fn=transcriptome)

#genmodpos=23510
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'position.pdf')

#with open("/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/corresp_tx_genome_tables/UAAUA_23510.txt", "r") as a_file:
#       for line in a_file:
#               split_string = line.split("\t")
#               modpos=int(split_string[0])
#               tx=a_string = split_string[1].rstrip("\n")
#               print(tx)
#               #p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,9])
#               fig, ax = s.plot_position (ref_id=tx, pos=modpos)
#               pdf.savefig(fig[0])
#pdf.close()

#kmer='GCCAU'
#s=SampCompDB(rootdir+"/outSampComp.db", fasta_fn=transcriptome)

#genmodpos=26337
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')

#with open("/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/corresp_tx_genome_tables/GCCAU_26337.txt", "r") as a_file:
#       for line in a_file:
#               split_string = line.split("\t")
#               modpos=int(split_string[0])
#               tx=a_string = split_string[1].rstrip("\n")
#               print(tx)
#               p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,9])
#               pdf.savefig(p1[0])
#pdf.close()

#kmer='CUUGA'
#s=SampCompDB(rootdir+"/outSampComp.db", fasta_fn=transcriptome)

#genmodpos=28960
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')

#with open("/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/corresp_tx_genome_tables/CUUGA_28960.txt", "r") as a_file:
#       for line in a_file:
#               split_string = line.split("\t")
#               modpos=int(split_string[0])
#               tx=a_string = split_string[1].rstrip("\n")
#               print(tx)
#               p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,9])
#               pdf.savefig(p1[0])
#pdf.close()

###############################################################################################################################################################################################################################################################################
                                                                                        ###     PUS7KD VS IVT (C37)      ###
###############################################################################################################################################################################################################################################################################

### UGUAR
#rootdir = '/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD_IVT/'
#resdir = '/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD_IVT/plots/'
#tabledir = '/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/corresp_tx_genome_tables/'
#os.mkdir(resdir)


#kmer='UGUAA'
#genmodpos=51
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')

##with open(tabledir+kmer+"_"+str(genmodpos)+".txt", "r") as a_file:
#	for line in a_file:
#		split_string = line.split("\t")
#		modpos=int(split_string[0])
#		tx=a_string = split_string[1].rstrip("\n")
#		print(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db")
#		s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db", fasta_fn=transcriptome)
#		p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,11])
#		pdf.savefig(p1[0])
#pdf.close()
#
#kmer='UGUAA'
#genmodpos=29079
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')
#
#with open(tabledir+kmer+"_"+str(genmodpos)+".txt", "r") as a_file:
#        for line in a_file:
#                split_string = line.split("\t")
#                modpos=int(split_string[0])
#                tx=a_string = split_string[1].rstrip("\n")
#                print(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db")
#                s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db", fasta_fn=transcriptome)
#                p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,11])
#                pdf.savefig(p1[0])
#pdf.close()
#
#kmer='UGUAG'
#genmodpos=29649
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')
#
#with open(tabledir+kmer+"_"+str(genmodpos)+".txt", "r") as a_file:
#        for line in a_file:
#                split_string = line.split("\t")
#                modpos=int(split_string[0])
#                tx=a_string = split_string[1].rstrip("\n")
#                print(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db")
#                s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db", fasta_fn=transcriptome)
#                p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,11])
#                pdf.savefig(p1[0])
#pdf.close()
#
#kmer='UGUAA'
#genmodpos=29690
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')
#
#with open(tabledir+kmer+"_"+str(genmodpos)+".txt", "r") as a_file:
#        for line in a_file:
#                split_string = line.split("\t")
#                modpos=int(split_string[0])
#                tx=a_string = split_string[1].rstrip("\n")
#                print(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db")
#                s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db", fasta_fn=transcriptome)
#                p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,11])
#                pdf.savefig(p1[0])
#pdf.close()
#
#kmer='UGUAA'
#genmodpos=29808
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')
#
#with open(tabledir+kmer+"_"+str(genmodpos)+".txt", "r") as a_file:
#        for line in a_file:
#                split_string = line.split("\t")
#                modpos=int(split_string[0])
#                tx=a_string = split_string[1].rstrip("\n")
#                print(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db")
#                s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db", fasta_fn=transcriptome)
#                p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,11])
#                pdf.savefig(p1[0])
#pdf.close()




###############################################################################################################################################################################################################################################################################
                                                                                        ###     WT VS IVT (C37)      ###
###############################################################################################################################################################################################################################################################################

### UGUAR
#rootdir = '/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/WT_IVT/'
#resdir = '/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/WT_IVT/plots/'
#tabledir = '/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/corresp_tx_genome_tables/'
#os.mkdir(resdir)


#kmer='UGUAA'
#genmodpos=51
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')

#with open(tabledir+kmer+"_"+str(genmodpos)+".txt", "r") as a_file:
#       for line in a_file:
#               split_string = line.split("\t")
#               modpos=int(split_string[0])
#               tx=a_string = split_string[1].rstrip("\n")
#               print(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db")
#               s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db", fasta_fn=transcriptome)
#               p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,11])
#               pdf.savefig(p1[0])
#pdf.close()

#kmer='UGUAA'
#genmodpos=29079
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')

#with open(tabledir+kmer+"_"+str(genmodpos)+".txt", "r") as a_file:
#        for line in a_file:
#                split_string = line.split("\t")
#                modpos=int(split_string[0])
#                tx=a_string = split_string[1].rstrip("\n")
#                print(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db")
#                s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db", fasta_fn=transcriptome)
#                p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,11])
#                pdf.savefig(p1[0])
#pdf.close()

#kmer='UGUAG'
#genmodpos=29649
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')

#with open(tabledir+kmer+"_"+str(genmodpos)+".txt", "r") as a_file:
#        for line in a_file:
#                split_string = line.split("\t")
#                modpos=int(split_string[0])
#                tx=a_string = split_string[1].rstrip("\n")
#                print(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db")
#                s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db", fasta_fn=transcriptome)
#                p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,11])
#                pdf.savefig(p1[0])
#pdf.close()

#kmer='UGUAA'
#genmodpos=29690
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')

#with open(tabledir+kmer+"_"+str(genmodpos)+".txt", "r") as a_file:
#        for line in a_file:
#                split_string = line.split("\t")
#                modpos=int(split_string[0])
#                tx=a_string = split_string[1].rstrip("\n")
#                print(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db")
#                s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db", fasta_fn=transcriptome)
#                p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,11])
#                pdf.savefig(p1[0])
#pdf.close()

#kmer='UGUAA'
#genmodpos=29808
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'_'+str(genmodpos)+'.pdf')

#with open(tabledir+kmer+"_"+str(genmodpos)+".txt", "r") as a_file:
#        for line in a_file:
#                split_string = line.split("\t")
#                modpos=int(split_string[0])
#                tx=a_string = split_string[1].rstrip("\n")
#                print(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db")
#                s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db", fasta_fn=transcriptome)
#                p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,11])
#                pdf.savefig(p1[0])
#pdf.close()


###############################################################################################################################################################################################################################################################################
                                                                                        ###     WT VS Pus7KD (C37) 	HUMAN TRANSCRIPTS      ###
###############################################################################################################################################################################################################################################################################

### UGUAR
#rootdir = '/hpcnfs/scratch/TSSM/tleonardi/synthetic_oligos/nanocompore_pipeline/results/nanocompore_oligo2'
#resdir = '/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_human_transcriptome/NANOCOMPORE/oligo2_plots/'
#os.mkdir(resdir)

#kmer='UGUAG'
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+'.pdf')
#s=SampCompDB(rootdir+"/out_SampComp.db", fasta_fn="/hpcnfs/scratch/TSSM/tleonardi/synthetic_oligos/nanocompore_pipeline/results/references/reference_transcriptome.fa")
#p1 = s.plot_signal(ref_id="control", start=60, end=70, plot_style='seaborn-whitegrid', figsize=[18,11])
#pdf.savefig(p1[0])
#pdf.close()


###############################################################################################################################################################################################################################################################################
                                                                                        ###       KIM VS IVT      ###
###############################################################################################################################################################################################################################################################################


#rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/nanocompore/comparison/kim/'
#resdir = '/hpcnfs/scratch/temporary/camilla_TL/5p_analysis/kim_vs_IVT/'
#tabledir = '/hpcnfs/scratch/temporary/camilla_TL/5p_analysis/kim_vs_IVT/corresp_tables/'
#os.mkdir(resdir)

#kmer='ACUUU'
#genmodpos=38
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+"_"+str(genmodpos)+'.pdf')

#with open(tabledir+kmer+"_"+str(genmodpos)+".txt", "r") as a_file:
#       for line in a_file:
#               split_string = line.split("\t")
#               modpos=int(split_string[0])
#               tx=a_string = split_string[1].rstrip("\n")
#               print(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/out_SampComp.db")
#               s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/out_SampComp.db", fasta_fn=transcriptome)
#               p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,14])
#               pdf.savefig(p1[0])
#pdf.close()




#rootdir = '/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD_WT/'
#resdir = '/hpcnfs/scratch/temporary/camilla_TL/5p_analysis/wt_vs_pus7kd/'
#tabledir = '/hpcnfs/scratch/temporary/camilla_TL/5p_analysis/kim_vs_IVT/corresp_tables/'
#os.mkdir(resdir)

#s=SampCompDB(rootdir+"/outSampComp.db", fasta_fn=transcriptome)
#kmer='ACUUU'
#genmodpos=38
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+kmer+"_"+str(genmodpos)+'.pdf')

#with open(tabledir+kmer+"_"+str(genmodpos)+".txt", "r") as a_file:
#       for line in a_file:
#               split_string = line.split("\t")
#               modpos=int(split_string[0])
#               tx=a_string = split_string[1].rstrip("\n")
#               print(tx)
#               p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,14])
#               pdf.savefig(p1[0])
#pdf.close()


