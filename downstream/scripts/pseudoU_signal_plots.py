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

### SIGNALS PLOTS 

transcriptome="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
#fasta_sequences = SeqIO.parse(open(transcriptome),'fasta')
#for fasta in fasta_sequences:
#	print(fasta.id)

## create a pdf for each sample (every sample has n plots = n transcripts)

###############################################################################################################################################################################################################################################################################
											###	High-confidence U sites 	###
###############################################################################################################################################################################################################################################################################
#rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/vero/NANOCOMPORE/sampcomp/'
#resdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/vero/NANOCOMPORE/sampcomp/plots/'
#tabledir="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/downstream/results_allfiles_LOR05_pval001/shared/"
#os.mkdir(resdir)
#for file in os.listdir(tabledir):
#	if file.endswith('txt'):
#		filename = os.path.join(tabledir, file)
#		genmodpos=file.split('_')[0]
#		kmer=file.split('_')[1].split('.')[0]
#		print(genmodpos+kmer)
#		pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos)+'_'+kmer+'.pdf')
#
#		with open(filename, "r") as a_file:
#			for line in a_file:
#				split_string = line.split("\t")
#				modpos=int(split_string[0])
#				tx=a_string = split_string[1].rstrip("\n")
#				s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db", fasta_fn=transcriptome)
#				p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,9])
#				pdf.savefig(p1[0])
#		pdf.close()



#########	BURROWS SITES in all cell lines

cell_line_list = [ 'caco2', 'calu3']

for cell_line in cell_line_list:
	print(cell_line)
	rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/'+cell_line+'/NANOCOMPORE/sampcomp/'
	resdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/'+cell_line+'/NANOCOMPORE/sampcomp/plots/Burrows/'
	tabledir='/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/downstream/results_allfiles_LOR05_pval001/corresp_tables_Burrows/'+cell_line
	if not os.path.exists(resdir):
		os.mkdir(resdir)
	for file in os.listdir(tabledir):
		if file.endswith('txt'):
			filename = os.path.join(tabledir, file)
			genmodpos=file.split('_')[0]
			kmer=file.split('_')[1].split('.')[0]
			print(genmodpos+kmer)
			pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos)+'_'+kmer+'.pdf')

			with open(filename, "r") as a_file:
				for line in a_file:
					split_string = line.split("\t")
					modpos=int(split_string[0])
					tx=a_string = split_string[1].rstrip("\n")
					s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db", fasta_fn=transcriptome)
					p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,9])
					pdf.savefig(p1[0])
			pdf.close()
	



