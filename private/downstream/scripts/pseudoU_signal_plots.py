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
#tabledir="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/downstream/results_allfiles_LOR05_pval001/shared/corresp_tables/"
#os.mkdir(resdir)
#cell_line_list = ['calu3']

#for cell_line in cell_line_list:
#	print(cell_line)
#	rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/nanocompore/comparison/sraf_calu3/'
#	resdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/'+cell_line+'/NANOCOMPORE/sampcomp/peakcalled_plots/'
#	tabledir='/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/downstream/results_allfiles_LOR05_pval001/shared/corresp_tables/'
#	if not os.path.exists(resdir):
#		os.mkdir(resdir)

#	for file in os.listdir(tabledir):
#		if file.endswith('txt'):
#			filename = os.path.join(tabledir, file)
#			genmodpos=file.split('_')[0]
#			kmer=file.split('_')[1].split('.')[0]
#			print(genmodpos+kmer)
#			pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos)+'_'+kmer+'.pdf')

#			with open(filename, "r") as a_file:
#				for line in a_file:
#					split_string = line.split("\t")
#					modpos=int(split_string[0])
#					tx=a_string = split_string[1].rstrip("\n")
#					s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/out_SampComp.db", fasta_fn=transcriptome)
#					p1 = s.plot_signal(ref_id=tx, start=(modpos-15), end=(modpos+15), plot_style='seaborn-whitegrid', figsize=[18,9])
#					pdf.savefig(p1[0])
#			pdf.close()



############################################	oligo2 Nanocompore paper ##############################################

#rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/data/nanocompore_paper_oligos/'
#resdir = '/hpcnfs/scratch/TSSM/cugolini/cov/data/nanocompore_paper_oligos/'
#reference=rootdir+"/ref.fa"
#modpos=65
#tx='control'
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(modpos)+'_'+'.pdf')
#s=SampCompDB(rootdir+"/oligo2SampComp.db", fasta_fn=reference)
#p1 = s.plot_signal(ref_id=tx, start=(modpos-15), end=(modpos+15), plot_style='seaborn-whitegrid', figsize=[18,9])
#pdf.savefig(p1[0])
#pdf.close()

###############################################################################################################################################################################################################################################################################	Paper sites ########################################################################################################################################################################################################################################################################################################################################

###     SPIKE

cell_line_list = ['caco2','vero']
for cell_line in cell_line_list:
	print(cell_line)
	rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/'+cell_line+'/NANOCOMPORE/sampcomp/'
	resdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/'+cell_line+'/NANOCOMPORE/sampcomp/paper_plots/'
	if not os.path.exists(resdir):
		os.mkdir(resdir)

	####### SPIKE
	#tx='3f2f5640-da4a-4237-952f-cd4a41e84234_315_NC_045512v2_11-29871'
	#tx1='3f2f5640-da4a-4237-952f-cd4a41e84234|315::NC_045512v2:11-29871'

	#genmodpos_start=23622
	#genmodpos_end=23638
	#modpos_start=2120
	#modpos_end=2141

	#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	#s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	#p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[18,9])
	#pdf.savefig(p1[0])
	#pdf.close()

	###### VME1
	#tx='fc4d6218-3c79-466d-84e4-705bee2dacef_13813_NC_045512v2_15-29873'
	#tx1='fc4d6218-3c79-466d-84e4-705bee2dacef|13813::NC_045512v2:15-29873'

	genmodpos_start=37
	genmodpos_end=54
	#modpos_start=22
	#modpos_end=46

	pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	#s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	#p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[18,9])
	#pdf.savefig(p1[0])

	tx='ade569a0-3419-462e-b6c5-77d6f75a8ed5_263_NC_045512v2_11-29871'
	tx1='ade569a0-3419-462e-b6c5-77d6f75a8ed5|263::NC_045512v2:11-29871'

	modpos_start=21
	modpos_end=43

	s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[18,9])
	pdf.savefig(p1[0])

	#tx='4167c2ea-24da-4951-a5e2-58927d55509a_143_NC_045512v2_14-29871'
	#tx1='4167c2ea-24da-4951-a5e2-58927d55509a|143::NC_045512v2:14-29871'

	#modpos_start=23
	#modpos_end=47

	#s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	#p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[18,9])
	#pdf.savefig(p1[0])
	pdf.close()


	#######	AP3A
	#tx='e94116d8-72eb-4252-819c-ad06fd9870d4_4145_NC_045512v2_14-29871'
	#tx1='e94116d8-72eb-4252-819c-ad06fd9870d4|4145::NC_045512v2:14-29871'

	#genmodpos_start=26391
	#genmodpos_end=26399
	#modpos_start=1057
	#modpos_end=1070

	#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	#s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	#p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[18,9])
	#pdf.savefig(p1[0])	
	#pdf.close()


	#genmodpos_start=29652
	#genmodpos_end=29660
	#modpos_start=4318
	#modpos_end=4331

	#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	#p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[18,9])
	#pdf.savefig(p1[0])
	#pdf.close()

        ####### NS6
	#tx='e94116d8-72eb-4252-819c-ad06fd9870d4_4145_NC_045512v2_14-29871'
	#tx1='e94116d8-72eb-4252-819c-ad06fd9870d4|4145::NC_045512v2:14-29871'

	#genmodpos_start=47
	#genmodpos_end=55
	#modpos_start=28
	#modpos_end=42

	#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	#s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	#p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[18,9])
	#pdf.savefig(p1[0])
	#pdf.close()

        ####### ORF10
	#tx='efad7b96-ac2e-4ce1-9b83-863ffdb18eac_116_NC_045512v2_11-29873'
	#tx1='efad7b96-ac2e-4ce1-9b83-863ffdb18eac|116::NC_045512v2:11-29873'

	#genmodpos_start=41
	#genmodpos_end=54
	#modpos_start=27
	#modpos_end=44

	#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	#s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	#p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[18,9])
	#pdf.savefig(p1[0])
	#pdf.close()

	#genmodpos_start=29640
	#genmodpos_end=29678
	#modpos_start=134
	#modpos_end=175

	#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	#s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	#p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[25,9])
	#pdf.savefig(p1[0])
	#pdf.close()

	#genmodpos_start=29690
	#genmodpos_end=29698
	#modpos_start=183
	#modpos_end=195

	#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	#s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	#p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[25,9])
	#pdf.savefig(p1[0])
	#pdf.close()



########################à	peakcalled sites on multiple isoforms

###genPos=29845
#vero
#rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/vero/NANOCOMPORE/sampcomp/'
#resdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/vero/NANOCOMPORE/sampcomp/peakcalled_plots/'
#genmodpos=29845
#modpos=1641
#kmer='UUUAA'
#tx='61a628be-eac5-4042-9c08-9674a992a139_1817_NC_045512v2_14-29873'
#tx1='61a628be-eac5-4042-9c08-9674a992a139|1817::NC_045512v2:14-29873'
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos)+'_'+kmer+'.pdf')
#s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
#p1 = s.plot_signal(ref_id=tx1, start=(modpos-15), end=(modpos+15), plot_style='seaborn-whitegrid', figsize=[18,9])
#pdf.savefig(p1[0])
#pdf.close()

#caco2
#rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/caco2/NANOCOMPORE/sampcomp/'
#resdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/caco2/NANOCOMPORE/sampcomp/peakcalled_plots/'
#genmodpos=29845
#modpos=1641
#kmer='UUUAA'
#tx='61a628be-eac5-4042-9c08-9674a992a139_1817_NC_045512v2_14-29873'
#tx1='61a628be-eac5-4042-9c08-9674a992a139|1817::NC_045512v2:14-29873'
#pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos)+'_'+kmer+'.pdf')
#s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
#p1 = s.plot_signal(ref_id=tx1, start=(modpos-15), end=(modpos+15), plot_style='seaborn-whitegrid', figsize=[18,9])
#pdf.savefig(p1[0])
#pdf.close()
'''
#calu3
rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/nanocompore/comparison/sraf_calu3/'
resdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/calu3/NANOCOMPORE/sampcomp/peakcalled_plots/'
genmodpos=29845
modpos=1641
kmer='UUUAA'
tx='61a628be-eac5-4042-9c08-9674a992a139_1817_NC_045512v2_14-29873'
tx1='61a628be-eac5-4042-9c08-9674a992a139|1817::NC_045512v2:14-29873'
pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos)+'_'+kmer+'.pdf')
s=SampCompDB(rootdir+tx+"/out_SampComp.db", fasta_fn=transcriptome)
p1 = s.plot_signal(ref_id=tx1, start=(modpos-15), end=(modpos+15), plot_style='seaborn-whitegrid', figsize=[18,9])
pdf.savefig(p1[0])
pdf.close()

###genPos=51
#vero
rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/vero/NANOCOMPORE/sampcomp/'
resdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/vero/NANOCOMPORE/sampcomp/peakcalled_plots/'
genmodpos=48
modpos=34
kmer='UCUUG'
tx='5820681b-6191-45dd-959c-96a8097aafdd_2255_NC_045512v2_14-29871'
tx1='5820681b-6191-45dd-959c-96a8097aafdd|2255::NC_045512v2:14-29871'
pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos)+'_'+kmer+'.pdf')
s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
p1 = s.plot_signal(ref_id=tx1, start=(modpos-15), end=(modpos+15), plot_style='seaborn-whitegrid', figsize=[18,9])
pdf.savefig(p1[0])
pdf.close()

#caco2
rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/caco2/NANOCOMPORE/sampcomp/'
resdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/caco2/NANOCOMPORE/sampcomp/peakcalled_plots/'
genmodpos=48
modpos=34
kmer='UCUUG'
tx='5820681b-6191-45dd-959c-96a8097aafdd_2255_NC_045512v2_14-29871'
tx1='5820681b-6191-45dd-959c-96a8097aafdd|2255::NC_045512v2:14-29871'
pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos)+'_'+kmer+'.pdf')
s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
p1 = s.plot_signal(ref_id=tx1, start=(modpos-15), end=(modpos+15), plot_style='seaborn-whitegrid', figsize=[18,9])
pdf.savefig(p1[0])
pdf.close()


#calu3
rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/nanocompore/comparison/sraf_calu3/'
resdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/calu3/NANOCOMPORE/sampcomp/peakcalled_plots/'
genmodpos=48
modpos=34
kmer='UCUUG'
tx='5820681b-6191-45dd-959c-96a8097aafdd_2255_NC_045512v2_14-29871'
tx1='5820681b-6191-45dd-959c-96a8097aafdd|2255::NC_045512v2:14-29871'
pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos)+'_'+kmer+'.pdf')
s=SampCompDB(rootdir+tx+"/out_SampComp.db", fasta_fn=transcriptome)
p1 = s.plot_signal(ref_id=tx1, start=(modpos-15), end=(modpos+15), plot_style='seaborn-whitegrid', figsize=[18,9])
pdf.savefig(p1[0])
pdf.close()
'''



'''
###	SPIKE
rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/vero/NANOCOMPORE/sampcomp/'
resdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/vero/NANOCOMPORE/sampcomp/peakcalled_plots/'
tx='3f2f5640-da4a-4237-952f-cd4a41e84234_315_NC_045512v2_11-29871'
tx1='3f2f5640-da4a-4237-952f-cd4a41e84234|315::NC_045512v2:11-29871'

genmodpos_start=21852
genmodpos_end=21862
modpos_start=354
modpos_end=362

pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'.pdf')
s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
p1 = s.plot_signal(ref_id=tx1, start=(modpos_start-10), end=(modpos_end+10), plot_style='seaborn-whitegrid', figsize=[18,9])
pdf.savefig(p1[0])
pdf.close()


genmodpos_start=22169
genmodpos_end=22179
modpos_start=671
modpos_end=679

pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'.pdf')
s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
p1 = s.plot_signal(ref_id=tx1, start=(modpos_start-10), end=(modpos_end+10), plot_style='seaborn-whitegrid', figsize=[18,9])
pdf.savefig(p1[0])
pdf.close()


genmodpos_start=23620
genmodpos_end=23632
modpos_start=2122
modpos_end=2132

pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'.pdf')
s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
p1 = s.plot_signal(ref_id=tx1, start=(modpos_start-10), end=(modpos_end+10), plot_style='seaborn-whitegrid', figsize=[18,9])
pdf.savefig(p1[0])
pdf.close()


genmodpos_start=29566
genmodpos_end=29578
modpos_start=8068
modpos_end=8078

pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'.pdf')
s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
p1 = s.plot_signal(ref_id=tx1, start=(modpos_start-10), end=(modpos_end+10), plot_style='seaborn-whitegrid', figsize=[18,9])
pdf.savefig(p1[0])
pdf.close()



###     VME1
#vero
tx='fc4d6218-3c79-466d-84e4-705bee2dacef_13813_NC_045512v2_15-29873'
tx1='fc4d6218-3c79-466d-84e4-705bee2dacef|13813::NC_045512v2:15-29873'

genmodpos_start=28013
genmodpos_end=28025
modpos_start=1593
modpos_end=1603

pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+tx+'.pdf')
s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
p1 = s.plot_signal(ref_id=tx1, start=(modpos_start-10), end=(modpos_end+10), plot_style='seaborn-whitegrid', figsize=[18,9])
pdf.savefig(p1[0])
pdf.close()


genmodpos_start=29842
genmodpos_end=29852
modpos_start=3422
modpos_end=3430

pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+tx+'.pdf')
s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
p1 = s.plot_signal(ref_id=tx1, start=(modpos_start-10), end=(modpos_end+10), plot_style='seaborn-whitegrid', figsize=[18,9])
pdf.savefig(p1[0])
pdf.close()


tx='ade569a0-3419-462e-b6c5-77d6f75a8ed5_263_NC_045512v2_11-29871'
tx1='ade569a0-3419-462e-b6c5-77d6f75a8ed5|263::NC_045512v2:11-29871'

genmodpos_start=28013
genmodpos_end=28025
modpos_start=1596
modpos_end=1606

pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+tx+'.pdf')
s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
p1 = s.plot_signal(ref_id=tx1, start=(modpos_start-10), end=(modpos_end+10), plot_style='seaborn-whitegrid', figsize=[18,9])
pdf.savefig(p1[0])
pdf.close()


genmodpos_start=29842
genmodpos_end=29852
modpos_start=3425
modpos_end=3433

pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+tx+'.pdf')
s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
p1 = s.plot_signal(ref_id=tx1, start=(modpos_start-10), end=(modpos_end+10), plot_style='seaborn-whitegrid', figsize=[18,9])
pdf.savefig(p1[0])
pdf.close()

###     ORF9D
rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/vero/NANOCOMPORE/sampcomp/'
resdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/vero/NANOCOMPORE/sampcomp/peakcalled_plots/'
tx='de81ef19-655d-4ced-a9cc-cb8384001058_107_NC_045512v2_11-29874'
tx1='de81ef19-655d-4ced-a9cc-cb8384001058|107::NC_045512v2:11-29874'

genmodpos_start=29413
genmodpos_end=29424
modpos_start=334
modpos_end=344

pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_fleming.pdf')
s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
p1 = s.plot_signal(ref_id=tx1, start=(modpos_start-10), end=(modpos_end+10), plot_style='seaborn-whitegrid', figsize=[18,9])
pdf.savefig(p1[0])
pdf.close()
'''

#########	BURROWS SITES in all cell lines

#cell_line_list = [ 'caco2', 'calu3']

#for cell_line in cell_line_list:
#	print(cell_line)
#	rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/'+cell_line+'/NANOCOMPORE/sampcomp/'
#	resdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/'+cell_line+'/NANOCOMPORE/sampcomp/plots/Burrows/'
#	tabledir='/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/downstream/results_allfiles_LOR05_pval001/corresp_tables_Burrows/'+cell_line
#	if not os.path.exists(resdir):
#		os.mkdir(resdir)
#	for file in os.listdir(tabledir):
#		if file.endswith('txt'):
#			filename = os.path.join(tabledir, file)
#			genmodpos=file.split('_')[0]
#			kmer=file.split('_')[1].split('.')[0]
#			print(genmodpos+kmer)
#			pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos)+'_'+kmer+'.pdf')
#
#			with open(filename, "r") as a_file:
#				for line in a_file:
#					split_string = line.split("\t")
#					modpos=int(split_string[0])
#					tx=a_string = split_string[1].rstrip("\n")
#					s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/outSampComp.db", fasta_fn=transcriptome)
#					p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,9])
#					pdf.savefig(p1[0])
#			pdf.close()
	


#rootdir = '/hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/nanocompore/comparison/kim/'
#resdir = '/hpcnfs/scratch/temporary/camilla_TL/ciao'
#tabledir='/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/downstream/results_allfiles_LOR05_pval001/corresp_tables_Fleming/vero/'
#if not os.path.exists(resdir):
#	os.mkdir(resdir)
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
#				s=SampCompDB(rootdir+re.sub(":","_",re.sub("::","_",re.sub("\\|","_",tx)))+"/out_SampComp.db", fasta_fn=transcriptome)
#				p1 = s.plot_signal(ref_id=tx, start=(modpos-10), end=(modpos+10), plot_style='seaborn-whitegrid', figsize=[18,9])
#				pdf.savefig(p1[0])
#				pdf.close()


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


