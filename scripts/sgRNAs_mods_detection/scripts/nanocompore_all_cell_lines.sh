#!/bin/bash
set -e -o pipefail

# load variables from general configuration file
if [ -z "$path" ]
then
      CURR_DIR=$1
else
      CURR_DIR=$path
fi

CONFIG=$(echo $CURR_DIR | rev | cut -d'/' -f3- |rev)                                    # obtain configuration file directory
source $CONFIG/general/config.sh

# load local configuration file
source $CURR_DIR/config.sh $CURR_DIR
# load images
source $CURR_DIR/images.sh $CURR_DIR

# take sample file for basecalling version
SAMPLE_FILE="${FILES}/${condition}_samples_${BASECALLING}.txt"
FASTA="$BASEDIR/analysis/fasta/$BASECALLING"
ALIGNMENTS="$BASEDIR/analysis/alignments/$BASECALLING/all_cell_lines"
WD="$BASEDIR/analysis/sgRNAs_mods_detection/$BASECALLING/all_cell_lines"


for selected_cell_line in CaCo2 CaLu3 VeroE6; do
        # copy all raw fast5 in a directory for each cell line
       mkdir -p $BASEDIR/analysis/fast5_total

       while IFS=$'\t' read sample fasta fast5 cell_line source; do
               cp $fast5/*.fast5 $BASEDIR/analysis/fast5_total
       done < <(grep $selected_cell_line $SAMPLE_FILE)
done

mkdir -p $WD
cat $FASTA/CaCo2/WT/*.fa $FASTA/CaLu3/WT/*.fa $FASTA/VeroE6/WT/*.fa > $FASTA/all_cell_lines/total.fa
$NANOCOMPORE f5c index -t $THREADS --iop $PARALLEL_JOBS -d $BASEDIR/analysis/fast5_total $FASTA/all_cell_lines/total.fa
$SINGC minimap2 -t $THREADS $TRANSCRIPTOME_PARAM $TRANSCRIPTOME_ASSEMBLY $FASTA/all_cell_lines/total.fa > $ALIGNMENTS/total.sam
$SINGC samtools view -h -F 2324 -Sb $ALIGNMENTS/total.sam | $SINGC samtools sort > $ALIGNMENTS/total.bam
$SINGC samtools index $ALIGNMENTS/total.bam
$NANOCOMPORE sh -c "f5c eventalign --rna --min-mapq 0 -t $THREADS -r $FASTA/all_cell_lines/total.fa -b $ALIGNMENTS/total.bam --g $TRANSCRIPTOME_ASSEMBLY --samples --print-read-names --scale-events --disable-cuda=yes --iop $PARALLEL_JOBS > $WD/eventalign/eventalign_total.txt"
$NANOCOMPORE sh -c "nanocompore eventalign_collapse -t $THREADS -i $WD/eventalign/eventalign_total.txt -o $WD/eventalign/collapse"
for filename in $WD/IVT/eventalign*; do
               base=${filename##*/eventalign_}
               $NANOCOMPORE nanocompore sampcomp --file_list1 $WD/eventalign/collapse/out_eventalign_collapse.tsv --file_list2 "$filename"/out_eventalign_collapse.tsv --label1 WT --label2 IVT --fasta $TRANSCRIPTOME_ASSEMBLY --outpath $WD/nanocompore/sampcomp/"$base" --overwrite --downsample_high_coverage 5000 --allow_warnings --min_coverage 30 --logit --nthreads $THREADS --bed $NANOCOMP_BED
done

