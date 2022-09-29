#!/bin/bash
set -e -o pipefail
#export LC_ALL=C

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
DATA="$BASEDIR/data"
SCRIPTDIR="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/public/gRNA_mods_detection/scripts"
FILES="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/public/gRNA_mods_detection/files"
WD="$BASEDIR/analysis/gRNA_mods"
IMG="$BASEDIR/img"
MOUNT_DIR="/hpcnfs/scratch"
REF_DATA="/hpcnfs/scratch/FN/camilla/nanopore/data"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
GENOME_PARAM="-k 8 -w 1 -ax splice -g 30000 -G 30000 -A1 -B2 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -un -p 0.7"
THREADS=10
TEMP_DIR="/hpcnfs/scratch/temporary/camilla_FN/gRNA_mods"

TRANSCRIPTOME_ASSEMBLY="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
NANOCOMP_FA="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/consensus_extracted.fa"
NANOCOMP_BED="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus_name_commas.bed"
nanocomp104="/hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_v1.0.4.sif"
#singularity pull tleonardi/nanocompore:v1.0.3
singpore="/hpcnfs/scratch/TSSM/cugolini/tools/porechop/porechop.simg"


# pull image
if [ ! -f "$IMG/nrceq_pipeline_latest.sif" ]; then
        cd $IMG
        singularity pull docker://cugolini/nrceq_pipeline:latest
fi

# singularity command
SINGC="singularity exec -B $MOUNT_DIR $IMG/nrceq_pipeline_latest.sif"


# select genomic RNAs from all the samples and merge them
mkdir -p $TEMP_DIR/fast5

#while IFS=$'\t' read -r sample sample_fasta sample_fast5 cell_line source; do

	#mkdir -p $WD/$cell_line/fastq
	#mkdir -p $WD/$cell_line/alignments_to_genome
	#if [ -d $sample_fasta ]; then
	#	cat $sample_fasta/*.fastq > $WD/$cell_line/fastq/"$sample".fastq
	#	awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' $WD/$cell_line/fastq/"$sample".fastq > $WD/$cell_line/fastq/"$sample".fa
	#	rm $WD/$cell_line/fastq/"$sample".fastq
	#elif [ -f $sample_fasta ]; then
	#	cp $sample_fasta $WD/$cell_line/fastq/"$sample".fa
	#fi
	#cat $WD/$cell_line/fastq/"$sample".fa >> $TEMP_DIR/all_samples.fa
	#$SINGC minimap2 $GENOME_PARAM -t $THREADS $GENOME_FA $WD/$cell_line/fastq/"$sample".fa* > $WD/$cell_line/alignments_to_genome/"$sample".sam
	#$SINGC samtools view -F 2068 -Sb $WD/$cell_line/alignments_to_genome/"$sample".sam > $WD/$cell_line/alignments_to_genome/"$sample".bam
	#$SINGC samtools sort $WD/$cell_line/alignments_to_genome/"$sample".bam > $WD/$cell_line/alignments_to_genome/"$sample"_sorted.bam
	#rm $WD/$cell_line/alignments_to_genome/"$sample".sam
	#rm $WD/$cell_line/alignments_to_genome/"$sample".bam
	#$SINGC samtools index $WD/$cell_line/alignments_to_genome/"$sample"_sorted.bam
	#$SINGC bedtools bamtobed -bed12 -i $WD/$cell_line/alignments_to_genome/"$sample"_sorted.bam > $WD/$cell_line/alignments_to_genome/"$sample"_sorted.bed
	#awk '{if($2<=45 && $3>=29850 && $10==1) print $4}' $WD/$cell_line/alignments_to_genome/"$sample"_sorted.bed > $WD/$cell_line/alignments_to_genome/"$sample"_gRNAs.txt
	#$SINGC python3 $SCRIPTDIR/extract_reads_bam.py $WD/$cell_line/alignments_to_genome/"$sample"_sorted.bam $WD/$cell_line/alignments_to_genome/"$sample"_gRNAs.bam $WD/$cell_line/alignments_to_genome/"$sample"_gRNAs.txt
	#rm $WD/$cell_line/alignments_to_genome/"$sample"_gRNAs.txt
	#cp $sample_fast5/*.fast5 $TEMP_DIR/fast5	# this line will be used for eventalign
#done < <(tail -n +2 $FILES/samples.txt)



mkdir -p $WD/all_cell_lines/alignments_to_genome
#bam_list=$(find $WD/*/alignments_to_genome/ -name "*_gRNAs.bam")
#$SINGC samtools merge $WD/all_cell_lines/alignments_to_genome/gRNAs.bam $bam_list

#$SINGC samtools view -h -F 2324 -Sb $WD/all_cell_lines/alignments_to_genome/gRNAs.bam | $SINGC samtools sort > $WD/all_cell_lines/alignments_to_genome/gRNAs_sorted.bam
#$SINGC samtools index $WD/all_cell_lines/alignments_to_genome/gRNAs_sorted.bam

# eventalign reads 

mkdir -p $WD/all_cell_lines/eventalign/collapse
$SINGC /f5c-v0.6/f5c_x86_64_linux index -t $THREADS --iop 5 -d $TEMP_DIR/fast5 $TEMP_DIR/all_samples.fa
#$SINGC f5c eventalign --rna --min-mapq 0 -t $THREADS -r $TEMP_DIR/all_samples.fa -b $WD/all_cell_lines/alignments_to_genome/gRNAs_sorted.bam --g $GENOME_FA --samples --print-read-names --scale-events --iop 5  | $SINGC NanopolishComp Eventalign_collapse -t $THREADS -o $WD/all_cell_lines/eventalign/collapse


#mkdir -p  $WD/caco2/NANOCOMPORE/eventalign/collapse/
#mkdir -p  $WD/caco2/NANOCOMPORE/sampcomp/
#FAST5_TOTAL="$WD/caco2/temporary_fast5/"
#PUBLIC_DATA="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi"



#for filename in /hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/eventalign/IVT/eventalign*; do
#       base=${filename##*/eventalign_}
#       singularity exec -B /hpcnfs/scratch/ $nanocomp104 nanocompore sampcomp --file_list1 $WD/caco2/NANOCOMPORE/eventalign/collapse/out_eventalign_collapse.tsv  --file_list2 "$filename"/out_eventalign_collapse.tsv  --label1 WT --label2 IVT --fasta $NANOCOMP_FA --outpath  $WD/caco2/NANOCOMPORE/sampcomp/"$base" --overwrite --downsample_high_coverage 5000 --allow_warnings --min_coverage 30 --logit --nthreads 10 --bed $NANOCOMP_BED
#done
