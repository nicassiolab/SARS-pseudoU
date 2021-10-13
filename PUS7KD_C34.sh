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
PUS7_KD="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210414_1407_X2_FAL77079_62af908f/S33250_CaCo2_C34_plus_doxy_Pus7KD"
WT="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210414_1407_X1_FAL82866_fb16777a/S33249_CaCo2_C34"
TRANSCRIPTOME_ASSEMBLY="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
NANOCOMP_FA="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/consensus_extracted.fa"
NANOCOMP_BED="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus_name_commas.bed"


#### 	SCRIPT FOR THE ANALYSIS OF PUS7 KD SAMPLES

# pycoQC

#source activate /hpcnfs/home/ieo5215/miniconda/envs/pycoQC
#mkdir -p $WD/PUS7_KD/PYCOQC
#pycoQC -f $PUS7_KD/sequencing_summary_FAL77079_12544e44.txt -o $WD/PYCOQC/sequencing_summary_PUS_KD.html
#pycoQC -f $WT/sequencing_summary_FAL82866_cf1b582e.txt -o $WD/PYCOQC/sequencing_summary_WT.html
#conda deactivate

# process fastq

mkdir -p $DATA/fastq_PUS7
#cat $PUS7_KD/fastq_pass/*.fastq > $DATA/fastq_PUS7/PUS7_KD.fastq
#cat $WT/fastq_pass/*.fastq > $DATA/fastq_PUS7/WT.fastq

# map to the assembly(EXTRACTION)

mkdir -p $WD/PUS7_KD/map_to_recap_assembly
#singularity exec -B /hpcnfs/scratch/ $singpore minimap2 -t 20 -ax map-ont -p 0 -N 10 $TRANSCRIPTOME_ASSEMBLY $DATA/fastq_PUS7/PUS7_KD.fastq > $WD/PUS7_KD/map_to_recap_assembly/PUS7_KD.sam
#singularity exec -B /hpcnfs/scratch/ $singpore minimap2 -t 20 -ax map-ont -p 0 -N 10 $TRANSCRIPTOME_ASSEMBLY $DATA/fastq_PUS7/WT.fastq > $WD/PUS7_KD/map_to_recap_assembly/WT.sam

###     CREATE ALIGNMENTS BACKUP, NANOCOUNT AND NANOCOMPORE FILE PREP

#mkdir -p $WD/PUS7_KD/map_to_recap_assembly/ALIGNMENTS_BACKUP
#mkdir -p $WD/PUS7_KD/map_to_recap_assembly/NANOCOUNT/counts
#mkdir -p $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome

#for filename in $WD/PUS7_KD/map_to_recap_assembly/*.sam; do
#        name=${filename##*/}
#        base=${name%.sam}
#        singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -Sb $filename > $WD/PUS7_KD/map_to_recap_assembly/ALIGNMENTS_BACKUP/"$base".bam
#        singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -F 2068 -Sb $filename > $WD/PUS7_KD/map_to_recap_assembly/NANOCOUNT/"$base"_filtered_np.bam
#        singularity exec -B /hpcnfs/scratch/ $singpore samtools sort $WD/PUS7_KD/map_to_recap_assembly/NANOCOUNT/"$base"_filtered_np.bam > $WD/PUS7_KD/map_to_recap_assembly/NANOCOUNT/"$base"_sorted_np.bam
#        singularity exec -B /hpcnfs/scratch/ $singpore samtools index $WD/PUS7_KD/map_to_recap_assembly/NANOCOUNT/"$base"_sorted_np.bam
#        rm $WD/PUS7_KD/map_to_recap_assembly/NANOCOUNT/"$base"_filtered_np.bam
#        singularity exec -B /hpcnfs/scratch/ $singtot NanoCount -3 10 -5 10 -p align_score -x -i $WD/PUS7_KD/map_to_recap_assembly/NANOCOUNT/"$base"_sorted_np.bam -o $WD/PUS7_KD/map_to_recap_assembly/NANOCOUNT/counts/"$base"_counts.tsv

#        singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -F 2324 -Sb $filename > $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/"$base"_filtered.bam
#        singularity exec -B /hpcnfs/scratch/ $singpore samtools sort $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/"$base"_filtered.bam > $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/"$base"_sorted.bam
#        singularity exec -B /hpcnfs/scratch/ $singpore samtools index  $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/"$base"_sorted.bam
#        rm $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/"$base"_filtered.bam
#        rm $filename
#done


### F5C

#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux index -d $PUS7_KD/fast5_pass $DATA/fastq_PUS7/PUS7_KD.fastq
#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux index -d $WT/fast5_pass $DATA/fastq_PUS7/WT.fastq

#mkdir -p $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign
#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux eventalign --rna --min-mapq 0 -t 20 -r $DATA/fastq_PUS7/PUS7_KD.fastq -b $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/PUS7_KD_sorted.bam --g $TRANSCRIPTOME_ASSEMBLY --samples --print-read-names --scale-events | singularity exec -B /hpcnfs/scratch/ /hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_8d6a70c.img  NanopolishComp Eventalign_collapse -o $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS7_KD/

#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux eventalign --rna --min-mapq 0 -t 20 -r $DATA/fastq_PUS7/WT.fastq -b $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/WT_sorted.bam --g $TRANSCRIPTOME_ASSEMBLY --samples --print-read-names --scale-events | singularity exec -B /hpcnfs/scratch/ /hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_8d6a70c.img  NanopolishComp Eventalign_collapse -o $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/WT

#source activate $ENVS/nanocompore_v100rc32/


### AGAINST IVT

#for filename in /hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/eventalign/IVT/eventalign*; do
#        base=${filename##*/eventalign_}
#	/hpcnfs/home/ieo5215/miniconda/envs/nanocompore_v100rc32/bin/nanocompore sampcomp --file_list1 $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS7_KD/out_eventalign_collapse.tsv --file_list2 "$filename"/out_eventalign_collapse.tsv  --label1 PUS7_KD --label2 IVT --fasta $NANOCOMP_FA --outpath $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD/"$base" --overwrite --downsample_high_coverage 5000 --allow_warnings --pvalue_thr 0.01 --min_coverage 30 --logit --nthreads 3 --bed $NANOCOMP_BED
#done


### WT AGAINST PUS7KD

#/hpcnfs/home/ieo5215/miniconda/envs/nanocompore_v100rc32/bin/nanocompore sampcomp --file_list1 $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/WT/out_eventalign_collapse.tsv --file_list2 $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS7_KD/out_eventalign_collapse.tsv  --label1 WT --label2 PUS7_KD --fasta $NANOCOMP_FA --outpath $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD_WT/ --overwrite --downsample_high_coverage 5000 --allow_warnings --pvalue_thr 0.01 --min_coverage 30 --logit --nthreads 3 --bed $NANOCOMP_BED

### reduce nanocompore filtering 
#singularity run -B /hpcnfs/scratch/ docker://tleonardi/nanocompore:v1.0.4 nanocompore sampcomp --file_list1 $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/WT/out_eventalign_collapse.tsv --file_list2 $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS7_KD/out_eventalign_collapse.tsv  --label1 WT --label2 PUS7_KD --fasta $NANOCOMP_FA --outpath $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD_WT_kmers_freq05/ --overwrite --max_invalid_kmers_freq 0.5  --downsample_high_coverage 5000 --allow_warnings --pvalue_thr 0.01 --min_coverage 30 --logit --nthreads 3 --bed $NANOCOMP_BED


### upsampling of KD (KD reads are much less than WT)
UPSAMPLED="$WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS7_KD_IVT_upsampled/"
#cp -r  $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS7_KD/ $UPSAMPLED
#tail -n +2 $UPSAMPLED/out_eventalign_collapse.tsv.idx > $UPSAMPLED/nohead
#for i in {1..6};do cat $UPSAMPLED/nohead >> $UPSAMPLED/temp; done
#cat $UPSAMPLED/out_eventalign_collapse.tsv.idx  $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/even/temp >  $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS_KD_upsampled/final
#rm $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS_KD_upsampled/out_eventalign_collapse.tsv.idx
#rm $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS_KD_upsampled/nohead
#rm $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS_KD_upsampled/temp
#mv $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS_KD_upsampled/final $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS_KD_upsampled/out_eventalign_collapse.tsv.idx
#mv $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS_KD_upsampled/ $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS7_KD_upsampled/


#singularity run -B /hpcnfs/scratch/ docker://tleonardi/nanocompore:v1.0.4 nanocompore sampcomp --file_list1 $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/WT/out_eventalign_collapse.tsv --file_list2 $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/PUS7_KD_upsampled/out_eventalign_collapse.tsv  --label1 WT --label2 PUS7_KD --fasta $NANOCOMP_FA --outpath $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/sampcomp/PUS7_KD_WT_upsampled/ --overwrite --max_invalid_kmers_freq 0.5  --downsample_high_coverage 5000 --allow_warnings --pvalue_thr 0.01 --min_coverage 30 --logit --nthreads 3 --bed $NANOCOMP_BED




############################################################################ map WT and KD to HUMAN TRANSCRIPTOME  ##########################################################################

mkdir -p $WD/PUS7_KD/map_to_human_transcriptome
mkdir -p $WD/PUS7_KD/map_to_human_transcriptome/ALIGNMENTS_BACKUP/
mkdir -p $WD/PUS7_KD/map_to_human_transcriptome/NANOCOUNT/counts

singularity exec -B /hpcnfs/scratch/ $singpore minimap2 -ax map-ont -uf -p 0 -N 10 -t 12 -I100g $REF_DATA/transcriptome_fasta_primary_assembly.fa $DATA/fastq_PUS7/WT.fastq > $WD/PUS7_KD/map_to_human_transcriptome/WT.sam
singularity exec -B /hpcnfs/scratch/ $singpore minimap2 -ax map-ont -uf -p 0 -N 10 -t 12 -I100g $REF_DATA/transcriptome_fasta_primary_assembly.fa $DATA/fastq_PUS7/PUS7_KD.fastq > $WD/PUS7_KD/map_to_human_transcriptome/PUS7_KD.sam
for filename in $WD/PUS7_KD/map_to_human_transcriptome/*.sam; do
        name=${filename##*/}
        base=${name%.sam}

        singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -Sb $filename > $WD/PUS7_KD/map_to_human_transcriptome/ALIGNMENTS_BACKUP/"$base".bam
        singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -F 2324 -b $WD/PUS7_KD/map_to_human_transcriptome/ALIGNMENTS_BACKUP/"$base".bam >  $WD/PUS7_KD/map_to_human_transcriptome/NANOCOMPORE/map_to_human_transcriptome/"$base".bam
        singularity exec -B /hpcnfs/scratch/ $singpore samtools sort $WD/PUS7_KD/map_to_human_transcriptome/NANOCOMPORE/map_to_human_transcriptome/"$base".bam > $WD/PUS7_KD/map_to_human_transcriptome/NANOCOMPORE/map_to_human_transcriptome/"$base"_sorted.bam
        singularity exec -B /hpcnfs/scratch/ $singpore samtools index $WD/PUS7_KD/map_to_human_transcriptome/NANOCOMPORE/map_to_human_transcriptome/"$base"_sorted.bam
        rm $WD/PUS7_KD/map_to_human_transcriptome/NANOCOMPORE/map_to_human_transcriptome/"$base".bam


        singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -F 2068 -Sb $filename > $WD/PUS7_KD/map_to_human_transcriptome/NANOCOUNT/"$base"_filtered_np.bam
        singularity exec -B /hpcnfs/scratch/ $singpore samtools sort $WD/PUS7_KD/map_to_human_transcriptome/NANOCOUNT/"$base"_filtered_np.bam > $WD/PUS7_KD/map_to_human_transcriptome/NANOCOUNT/"$base"_sorted_np.bam
        singularity exec -B /hpcnfs/scratch/ $singpore samtools index $WD/PUS7_KD/map_to_human_transcriptome/NANOCOUNT/"$base"_sorted_np.bam
        rm $WD/PUS7_KD/map_to_human_transcriptome/NANOCOUNT/"$base"_filtered_np.bam
        singularity exec -B /hpcnfs/scratch/ $singtot NanoCount -3 10 -5 10 -p align_score -x -i $WD/PUS7_KD/map_to_human_transcriptome/NANOCOUNT/"$base"_sorted_np.bam -o $WD/PUS7_KD/map_to_human_transcriptome/NANOCOUNT/counts/"$base"_counts.tsv
        rm $WD/PUS7_KD/map_to_human_transcriptome/"$base".sam

done


