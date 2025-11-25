#!/bin/python

import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
#from Bio import SeqIO
from nanocompore.SampCompDB import SampCompDB
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
	"--transcriptome_assembly",
	type=str
)
parser.add_argument(
	"--basecalling_version",
	type=str
)
parser.add_argument(
        "--base_directory",
        type=str
)
args = parser.parse_args()


transcriptome = args.transcriptome_assembly
basecalling = args.basecalling_version
BASEDIR = args.base_directory

### SIGNALS PLOTS 


cell_line_list = ['CaCo2','VeroE6']
for cell_line in cell_line_list:
	print(cell_line)
	rootdir = BASEDIR+'/analysis/sgRNAs_mods_detection/'+basecalling+"/"cell_line+'/nanocompore/sampcomp/'
	resdir = BASEDIR+'/results/sgRNAs_mods_detection/violin_plots/'+cell_line+'/'
	if not os.path.exists(resdir):
		os.mkdir(resdir)

	####### SPIKE
	tx='3f2f5640-da4a-4237-952f-cd4a41e84234_315_NC_045512v2_11-29871'
	tx1='3f2f5640-da4a-4237-952f-cd4a41e84234|315::NC_045512v2:11-29871'

	genmodpos_start=23622
	genmodpos_end=23638
	modpos_start=2120
	modpos_end=2141

	pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[18,9])
	pdf.savefig(p1[0])
	pdf.close()

	###### VME1
	tx='fc4d6218-3c79-466d-84e4-705bee2dacef_13813_NC_045512v2_15-29873'
	tx1='fc4d6218-3c79-466d-84e4-705bee2dacef|13813::NC_045512v2:15-29873'

	genmodpos_start=37
	genmodpos_end=54
	modpos_start=22
	modpos_end=46

	pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[18,9])
	pdf.savefig(p1[0])

	pdf.close()


	#######	AP3A
	tx='e94116d8-72eb-4252-819c-ad06fd9870d4_4145_NC_045512v2_14-29871'
	tx1='e94116d8-72eb-4252-819c-ad06fd9870d4|4145::NC_045512v2:14-29871'

	genmodpos_start=26391
	genmodpos_end=26399
	modpos_start=1057
	modpos_end=1070

	pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[18,9])
	pdf.savefig(p1[0])	
	pdf.close()


	genmodpos_start=29652
	genmodpos_end=29660
	modpos_start=4318
	modpos_end=4331

	pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[18,9])
	pdf.savefig(p1[0])
	pdf.close()

        ####### NS6
	tx='e94116d8-72eb-4252-819c-ad06fd9870d4_4145_NC_045512v2_14-29871'
	tx1='e94116d8-72eb-4252-819c-ad06fd9870d4|4145::NC_045512v2:14-29871'

	genmodpos_start=47
	genmodpos_end=55
	modpos_start=28
	modpos_end=42

	pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[18,9])
	pdf.savefig(p1[0])
	pdf.close()

        ####### ORF10
	tx='efad7b96-ac2e-4ce1-9b83-863ffdb18eac_116_NC_045512v2_11-29873'
	tx1='efad7b96-ac2e-4ce1-9b83-863ffdb18eac|116::NC_045512v2:11-29873'

	genmodpos_start=41
	genmodpos_end=54
	modpos_start=27
	modpos_end=44

	pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[18,9])
	pdf.savefig(p1[0])
	pdf.close()

	genmodpos_start=29640
	genmodpos_end=29678
	modpos_start=134
	modpos_end=175

	pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[25,9])
	pdf.savefig(p1[0])
	pdf.close()

	genmodpos_start=29690
	genmodpos_end=29698
	modpos_start=183
	modpos_end=195

	pdf = matplotlib.backends.backend_pdf.PdfPages(resdir+str(genmodpos_start)+'_'+str(genmodpos_end)+'_'+cell_line+'.pdf')
	s=SampCompDB(rootdir+tx+"/outSampComp.db", fasta_fn=transcriptome)
	p1 = s.plot_signal(ref_id=tx1, start=(modpos_start), end=(modpos_end), plot_style='seaborn-whitegrid', figsize=[25,9])
	pdf.savefig(p1[0])
	pdf.close()



