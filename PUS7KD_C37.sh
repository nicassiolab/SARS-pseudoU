#!/bin/bash
#set -e -o pipefail
#export LC_ALL=C

source /hpcnfs/home/ieo5215/miniconda/etc/profile.d/conda.sh

BASEDIR="/hpcnfs/scratch/FN/TL/cugolini/cov"
DATA="$BASEDIR/data"
WD="$BASEDIR/analysis"
REF_DATA="/hpcnfs/scratch/FN/camilla/nanopore/data"
ENVS="/hpcnfs/home/ieo5215/miniconda/envs"
singpore="/hpcnfs/scratch/TSSM/cugolini/tools/porechop/porechop.simg"
singfastp="/hpcnfs/scratch/TSSM/cugolini/tools/fastp/fastp.simg"
singtot="/hpcnfs/scratch/TSSM/cugolini/cov/img/recappable.simg"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
PUS7_KD="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210825_1506_X2_FAL77483_43bd693a/S35756_CaCo2_C37_plus_doxy_Pus7_Pus7L_KD"
WT="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210825_1506_X1_FAL77093_8c3d8e3b/S35755_CaCo2_C37"
TRANSCRIPTOME_ASSEMBLY="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
NANOCOMP_FA="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/consensus_extracted.fa"
NANOCOMP_BED="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus_name_commas.bed"

#### 	SCRIPT FOR THE ANALYSIS OF PUS7 KD SAMPLES

################## QC of the samples fastq
# pycoQC

#source activate /hpcnfs/home/ieo5215/miniconda/envs/pycoQC
#mkdir -p $WD/PUS7_KD_C37/PYCOQC
#pycoQC -f $PUS7_KD/sequencing_summary_FAL77483_f8ac5a2e.txt -o $WD/PUS7_KD_C37/PYCOQC/sequencing_summary_PUS_KD.html
#pycoQC -f $WT/sequencing_summary_FAL77093_5c5220a1.txt -o $WD/PUS7_KD_C37/PYCOQC/sequencing_summary_WT.html
#conda deactivate

# process fastq

#mkdir -p $DATA/fastq_PUS7_C37
#cat $PUS7_KD/fastq_pass/*.fastq > $DATA/fastq_PUS7_C37/PUS7_KD.fastq
#cat $WT/fastq_pass/*.fastq > $DATA/fastq_PUS7_C37/WT.fastq


############################################################################ map WT and KD to  VIRAL TRANSCRIPTOME (NRCEQ ASSEMBLY-EXTRACTION) ##########################################################################

# map to the assembly(EXTRACTION)
mkdir -p $WD/PUS7_KD_C37/map_to_recap_assembly
#singularity exec -B /hpcnfs/scratch/ $singpore minimap2 -t 20 -ax map-ont -p 0 -N 10 $TRANSCRIPTOME_ASSEMBLY $DATA/fastq_PUS7_C37/PUS7_KD.fastq > $WD/PUS7_KD_C37/map_to_recap_assembly/PUS7_KD.sam
#singularity exec -B /hpcnfs/scratch/ $singpore minimap2 -t 20 -ax map-ont -p 0 -N 10 $TRANSCRIPTOME_ASSEMBLY $DATA/fastq_PUS7_C37/WT.fastq > $WD/PUS7_KD_C37/map_to_recap_assembly/WT.sam

###     CREATE ALIGNMENTS BACKUP, NANOCOUNT AND NANOCOMPORE FILE PREP

mkdir -p $WD/PUS7_KD_C37/map_to_recap_assembly/ALIGNMENTS_BACKUP
mkdir -p $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOUNT/counts
mkdir -p $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome

#for filename in $WD/PUS7_KD_C37/map_to_recap_assembly/*.sam; do
#        name=${filename##*/}
#        base=${name%.sam}
#        singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -Sb $filename > $WD/PUS7_KD_C37/map_to_recap_assembly/ALIGNMENTS_BACKUP/"$base".bam
#        singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -F 2068 -Sb $filename > $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOUNT/"$base"_filtered_np.bam
#        singularity exec -B /hpcnfs/scratch/ $singpore samtools sort $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOUNT/"$base"_filtered_np.bam > $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOUNT/"$base"_sorted_np.bam
#        singularity exec -B /hpcnfs/scratch/ $singpore samtools index $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOUNT/"$base"_sorted_np.bam
#        rm $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOUNT/"$base"_filtered_np.bam
#        singularity exec -B /hpcnfs/scratch/ $singtot NanoCount -3 10 -5 10 -p align_score -x -i $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOUNT/"$base"_sorted_np.bam -o $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOUNT/counts/"$base"_counts.tsv

#        singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -F 2324 -Sb $filename > $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/"$base"_filtered.bam
#        singularity exec -B /hpcnfs/scratch/ $singpore samtools sort $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/"$base"_filtered.bam > $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/"$base"_sorted.bam
#        singularity exec -B /hpcnfs/scratch/ $singpore samtools index  $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/"$base"_sorted.bam
#        rm $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/"$base"_filtered.bam
#        rm $filename
#done


### F5C (eventalign)

#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux index -d $PUS7_KD/fast5_pass $DATA/fastq_PUS7_C37/PUS7_KD.fastq
#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux index -d $WT/fast5_pass $DATA/fastq_PUS7_C37/WT.fastq

mkdir -p $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/eventalign
#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux eventalign --rna --min-mapq 0 -t 20 -r $DATA/fastq_PUS7_C37/PUS7_KD.fastq -b $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/PUS7_KD_sorted.bam --g $TRANSCRIPTOME_ASSEMBLY --samples --print-read-names --scale-events | singularity exec -B /hpcnfs/scratch/ /hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_8d6a70c.img  NanopolishComp Eventalign_collapse -o $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS7_KD/

#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux eventalign --rna --min-mapq 0 -t 20 -r $DATA/fastq_PUS7_C37/WT.fastq -b $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/WT_sorted.bam --g $TRANSCRIPTOME_ASSEMBLY --samples --print-read-names --scale-events | singularity exec -B /hpcnfs/scratch/ /hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_8d6a70c.img  NanopolishComp Eventalign_collapse -o $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/eventalign/WT

#singularity pull docker://tleonardi/nanocompore:v1.0.4
#mv nanocompore_v1.0.4.sif /hpcnfs/scratch/TSSM/cugolini/cov/img/
nanocomp104="/hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_v1.0.4.sif"

#singularity exec -B /hpcnfs/scratch/ $nanocomp104 nanocompore sampcomp --file_list1 $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/eventalign/WT/out_eventalign_collapse.tsv --file_list2 $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS7_KD/out_eventalign_collapse.tsv  --label1 WT --label2 PUS7_KD --fasta $NANOCOMP_FA --outpath $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD_WT/ --overwrite --downsample_high_coverage 5000 --allow_warnings --pvalue_thr 0.01 --min_coverage 30 --logit --nthreads 3 --bed $NANOCOMP_BED

mkdir -p $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD_IVT/
mkdir -p $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/WT_IVT/

#for filename in /hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/eventalign/IVT/eventalign*; do
#        base=${filename##*/eventalign_}
#       singularity exec -B /hpcnfs/scratch/ $nanocomp104 nanocompore sampcomp --file_list1 $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS7_KD/out_eventalign_collapse.tsv --file_list2 "$filename"/out_eventalign_collapse.tsv  --label1 PUS7_KD --label2 IVT --fasta $NANOCOMP_FA --outpath $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD_IVT/"$base" --overwrite --downsample_high_coverage 5000 --allow_warnings --pvalue_thr 0.01 --min_coverage 30 --logit --nthreads 3 --bed $NANOCOMP_BED
#done

#for filename in /hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/eventalign/IVT/eventalign*; do
#        base=${filename##*/eventalign_}
#       singularity exec -B /hpcnfs/scratch/ $nanocomp104 nanocompore sampcomp --file_list1 $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/eventalign/WT/out_eventalign_collapse.tsv --file_list2 "$filename"/out_eventalign_collapse.tsv  --label1 WT --label2 IVT --fasta $NANOCOMP_FA --outpath $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/sampcomp/WT_IVT/"$base" --overwrite --downsample_high_coverage 5000 --allow_warnings --pvalue_thr 0.01 --min_coverage 30 --logit --nthreads 3 --bed $NANOCOMP_BED
#done




############################################################################ map WT and KD to  VIRAL GENOME ##########################################################################

mkdir -p $WD/PUS7_KD_C37/map_to_genome
#singularity exec -B /hpcnfs/scratch/ $singpore minimap2 -k 8 -w 1 -ax splice -g 30000 -G 30000 -A1 -B2 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -un -p 0.7 -t 8 $GENOME_FA $DATA/fastq_PUS7_C37/WT.fastq > $WD/PUS7_KD_C37/map_to_genome/WT.sam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -F 2324 -Sb $WD/PUS7_KD_C37/map_to_genome/WT.sam > $WD/PUS7_KD_C37/map_to_genome/WT.bam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools sort $WD/PUS7_KD_C37/map_to_genome/WT.bam > $WD/PUS7_KD_C37/map_to_genome/WT_sorted.bam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools index $WD/PUS7_KD_C37/map_to_genome/WT_sorted.bam
#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux eventalign --rna --min-mapq 0 -t 12 -r $DATA/fastq_PUS7_C37/WT.fastq -b $WD/PUS7_KD_C37/map_to_genome/WT_sorted.bam --g $GENOME_FA --samples --print-read-names --scale-events | singularity exec -B /hpcnfs/scratch/ /hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_8d6a70c.img  NanopolishComp Eventalign_collapse -o $WD/PUS7_KD_C37/map_to_genome/NANOCOMPORE/eventalign/WT

#singularity exec -B /hpcnfs/scratch/ $singpore minimap2 -k 8 -w 1 -ax splice -g 30000 -G 30000 -A1 -B2 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -un -p 0.7 -t 12 $GENOME_FA $DATA/fastq_PUS7_C37/PUS7_KD.fastq > $WD/PUS7_KD_C37/map_to_genome/PUS7_KD.sam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -F 2324 -Sb $WD/PUS7_KD_C37/map_to_genome/PUS7_KD.sam > $WD/PUS7_KD_C37/map_to_genome/PUS7_KD.bam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools sort $WD/PUS7_KD_C37/map_to_genome/PUS7_KD.bam > $WD/PUS7_KD_C37/map_to_genome/PUS7_KD_sorted.bam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools index $WD/PUS7_KD_C37/map_to_genome/PUS7_KD_sorted.bam
#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux eventalign --rna --min-mapq 0 -t 12 -r $DATA/fastq_PUS7_C37/PUS7_KD.fastq -b $WD/PUS7_KD_C37/map_to_genome/PUS7_KD_sorted.bam --g $GENOME_FA --samples --print-read-names --scale-events | singularity exec -B /hpcnfs/scratch/ /hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_8d6a70c.img  NanopolishComp Eventalign_collapse -o $WD/PUS7_KD_C37/map_to_genome/NANOCOMPORE/eventalign/PUS7_KD

#singularity exec -B /hpcnfs/scratch/ $nanocomp104 nanocompore sampcomp --file_list1 $WD/PUS7_KD_C37/map_to_genome/NANOCOMPORE/eventalign/WT/out_eventalign_collapse.tsv --file_list2 $WD/PUS7_KD_C37/map_to_genome/NANOCOMPORE/eventalign/PUS7_KD/out_eventalign_collapse.tsv  --label1 WT --label2 PUS7_KD --fasta $GENOME_FA --outpath $WD/PUS7_KD_C37/map_to_genome/NANOCOMPORE/sampcomp/PUS7_KD_WT/ --overwrite --downsample_high_coverage 5000 --allow_warnings --pvalue_thr 0.01 --min_coverage 30 --logit --nthreads 3

#IVT="/hpcnfs/scratch/TSSM/tleonardi/SARS-CoV-2-datasets/Kim_2020/fastq/IVT/ont_research.fa"

#singularity exec -B /hpcnfs/scratch/ $singpore minimap2 -k 8 -w 1 -ax splice -g 30000 -G 30000 -A1 -B2 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -un -p 0.7 -t 8 $GENOME_FA $IVT > $WD/PUS7_KD_C37/map_to_genome/IVT.sam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -F 2324 -Sb $WD/PUS7_KD_C37/map_to_genome/IVT.sam > $WD/PUS7_KD_C37/map_to_genome/IVT.bam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools sort $WD/PUS7_KD_C37/map_to_genome/IVT.bam > $WD/PUS7_KD_C37/map_to_genome/IVT_sorted.bam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools index $WD/PUS7_KD_C37/map_to_genome/IVT_sorted.bam
#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux eventalign --rna --min-mapq 0 -t 8 -r $IVT -b $WD/PUS7_KD_C37/map_to_genome/IVT_sorted.bam --g $GENOME_FA --samples --print-read-names --scale-events | singularity exec -B /hpcnfs/scratch/ /hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_8d6a70c.img  NanopolishComp Eventalign_collapse -o $WD/PUS7_KD_C37/map_to_genome/NANOCOMPORE/eventalign/IVT

#singularity exec -B /hpcnfs/scratch/ $nanocomp104 nanocompore sampcomp --file_list1 $WD/PUS7_KD_C37/map_to_genome/NANOCOMPORE/eventalign/WT/out_eventalign_collapse.tsv --file_list2 $WD/PUS7_KD_C37/map_to_genome/NANOCOMPORE/eventalign/IVT/out_eventalign_collapse.tsv  --label1 WT --label2 IVT --fasta $GENOME_FA --outpath $WD/PUS7_KD_C37/map_to_genome/NANOCOMPORE/sampcomp/WT_IVT/ --overwrite --downsample_high_coverage 5000 --allow_warnings --pvalue_thr 0.01 --min_coverage 30 --logit --nthreads 3

#singularity exec -B /hpcnfs/scratch/ $nanocomp104 nanocompore sampcomp --file_list1 $WD/PUS7_KD_C37/map_to_genome/NANOCOMPORE/eventalign/PUS7_KD/out_eventalign_collapse.tsv --file_list2 $WD/PUS7_KD_C37/map_to_genome/NANOCOMPORE/eventalign/IVT/out_eventalign_collapse.tsv  --label1 PUS7_KD --label2 IVT --fasta $GENOME_FA --outpath $WD/PUS7_KD_C37/map_to_genome/NANOCOMPORE/sampcomp/PUS7_KD_IVT/ --overwrite --downsample_high_coverage 5000 --allow_warnings --pvalue_thr 0.01 --min_coverage 30 --logit --nthreads 3




############################################################################ map WT and KD to HUMAN TRANSCRIPTOME  ##########################################################################

mkdir -p $WD/PUS7_KD_C37/map_to_human_transcriptome
mkdir -p $WD/PUS7_KD_C37/map_to_human_transcriptome/ALIGNMENTS_BACKUP/
mkdir -p $WD/PUS7_KD_C37/map_to_human_transcriptome/NANOCOUNT/counts

#singularity exec -B /hpcnfs/scratch/ $singpore minimap2 -ax map-ont -uf -p 0 -N 10 -t 20 -I100g $REF_DATA/transcriptome_fasta_primary_assembly.fa $DATA/fastq_PUS7_C37/WT.fastq > $WD/PUS7_KD_C37/map_to_human_transcriptome/WT.sam
#singularity exec -B /hpcnfs/scratch/ $singpore minimap2 -ax map-ont -uf -p 0 -N 10 -t 20 -I100g $REF_DATA/transcriptome_fasta_primary_assembly.fa $DATA/fastq_PUS7_C37/PUS7_KD.fastq > $WD/PUS7_KD_C37/map_to_human_transcriptome/PUS7_KD.sam
#for filename in $WD/PUS7_KD_C37/map_to_human_transcriptome/*.sam; do
        #name=${filename##*/}
        #base=${name%.sam}

        #singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -Sb $filename > $WD/PUS7_KD_C37/map_to_human_transcriptome/ALIGNMENTS_BACKUP/"$base".bam
	#singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -F 2324 -b $WD/PUS7_KD_C37/map_to_human_transcriptome/ALIGNMENTS_BACKUP/"$base".bam >  $WD/PUS7_KD_C37/map_to_human_transcriptome/NANOCOMPORE/map_to_human_transcriptome/"$base".bam
	#singularity exec -B /hpcnfs/scratch/ $singpore samtools sort $WD/PUS7_KD_C37/map_to_human_transcriptome/NANOCOMPORE/map_to_human_transcriptome/"$base".bam > $WD/PUS7_KD_C37/map_to_human_transcriptome/NANOCOMPORE/map_to_human_transcriptome/"$base"_sorted.bam 
	#singularity exec -B /hpcnfs/scratch/ $singpore samtools index $WD/PUS7_KD_C37/map_to_human_transcriptome/NANOCOMPORE/map_to_human_transcriptome/"$base"_sorted.bam
	#rm $WD/PUS7_KD_C37/map_to_human_transcriptome/NANOCOMPORE/map_to_human_transcriptome/"$base".bam
	

        #singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -F 2068 -Sb $filename > $WD/PUS7_KD_C37/map_to_human_transcriptome/NANOCOUNT/"$base"_filtered_np.bam
        #singularity exec -B /hpcnfs/scratch/ $singpore samtools sort $WD/PUS7_KD_C37/map_to_human_transcriptome/NANOCOUNT/"$base"_filtered_np.bam > $WD/PUS7_KD_C37/map_to_human_transcriptome/NANOCOUNT/"$base"_sorted_np.bam
        #singularity exec -B /hpcnfs/scratch/ $singpore samtools index $WD/PUS7_KD_C37/map_to_human_transcriptome/NANOCOUNT/"$base"_sorted_np.bam
        #rm $WD/PUS7_KD_C37/map_to_human_transcriptome/NANOCOUNT/"$base"_filtered_np.bam
        #singularity exec -B /hpcnfs/scratch/ $singtot NanoCount -3 10 -5 10 -p align_score -x -i $WD/PUS7_KD_C37/map_to_human_transcriptome/NANOCOUNT/"$base"_sorted_np.bam -o $WD/PUS7_KD_C37/map_to_human_transcriptome/NANOCOUNT/counts/"$base"_counts.tsv
	#rm $WD/PUS7_KD_C37/map_to_human_transcriptome/"$base".sam

#done


### F5C (eventalign)

#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux eventalign --rna --min-mapq 0 -t 20 -r $DATA/fastq_PUS7_C37/WT.fastq -b $WD/PUS7_KD_C37/map_to_human_transcriptome/WT_sorted.bam --g $REF_DATA/transcriptome_fasta_primary_assembly.fa --samples --print-read-names --scale-events | singularity exec -B /hpcnfs/scratch/ /hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_8d6a70c.img  NanopolishComp Eventalign_collapse -o $WD/PUS7_KD_C37/map_to_human_transcriptome/NANOCOMPORE/eventalign/WT
#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux eventalign --rna --min-mapq 0 -t 20 -r $DATA/fastq_PUS7_C37/PUS7_KD.fastq -b $WD/PUS7_KD_C37/map_to_human_transcriptome/PUS7_KD_sorted.bam --g $GENOME_FA --samples --print-read-names --scale-events | singularity exec -B /hpcnfs/scratch/ /hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_8d6a70c.img  NanopolishComp Eventalign_collapse -o $WD/PUS7_KD_C37/map_to_human_transcriptome/NANOCOMPORE/eventalign/PUS7_KD

