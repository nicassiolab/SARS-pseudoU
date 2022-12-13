#!/bin/bash
set -e -o pipefail

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
WD="$BASEDIR/analysis/sgRNAs_mods_detection"
FASTA="$BASEDIR/analysis/fasta"

# pull image for nanocompore
if [ ! -f "$IMG/nanocompore_v1.0.4.sif" ]; then
        cd $IMG
        singularity pull docker://tleonardi/nanocompore:v1.0.4
fi
# singularity command
NANOCOMPORE="singularity exec -B $MOUNT_DIR $IMG/nanocompore_v1.0.4.sif"

# indicate basecalling version (the first basecalling used in the analysis has mixed version therefore we just call it guppy_initial)
BASECALLING="guppy_v601"
WD="$WD/$BASECALLING"
FASTA="$FASTA/$BASECALLING"
ALIGNMENTS="$BASEDIR/analysis/alignments/$BASECALLING"

for condition in PUS7KD;do

	# take sample file for each condition and basecalling version
        SAMPLE_FILE="${FILES}/${condition}_samples_${BASECALLING}.txt"
	
	# read sample file 
        while IFS=$'\t' read sample fasta fast5 cell_line source; do
		SAMPLE="$sample"
                SAMPLE_FA="$fasta"
                SAMPLE_FAST5="$fast5"
                SAMPLE_CELL_LINE="$cell_line"
		mkdir -p $FASTA/$cell_line/"$condition"/


		$NANOCOMPORE f5c index --iop 5 -d $SAMPLE_FAST5 $FASTA/$cell_line/"$condition"/"$sample".fa
		mkdir -p $WD/$cell_line/"$condition"/alignments_to_assembly/eventalign/
		$NANOCOMPORE sh -c "f5c eventalign --rna --min-mapq 0 -t $THREADS -r $FASTA/$cell_line/"$condition"/"$sample".fa -b $ALIGNMENTS/$cell_line/"$condition"/alignments_to_assembly/"$sample"_nanocompore.bam --g $TRANSCRIPTOME_ASSEMBLY --samples --print-read-names --scale-events --disable-cuda=yes --iop 5 | nanocompore eventalign_collapse -o $WD/$cell_line/"$condition"/alignments_to_assembly/eventalign/collapse"



	done < <(tail -n +2 $SAMPLE_FILE)

done






# index and eventalign all the raw files from the same cell line together
SAMPLE_FILE="${FILES}/${condition_per_cell_line}_samples_${BASECALLING}.txt"

for selected_cell_line in CaCo2 VeroE6; do

        # copy all the raw fast5 files for each cell line into a directory
        mkdir -p $FAST5/$selected_cell_line

        while IFS=$'\t' read sample fasta fast5 cell_line source; do
                cp $fast5/*.fast5 $FAST5/$selected_cell_line
        done < <(grep $selected_cell_line $SAMPLE_FILE)

        # index and eventalign the files
        $NANOCOMPORE f5c index -t $THREADS -d $FAST5/$selected_cell_line $FASTA/per_cell_line/$condition_per_cell_line/"$selected_cell_line".fa
        $NANOCOMPORE f5c eventalign --rna --min-mapq 0 -t $THREADS -r $FASTA/per_cell_line/$condition_per_cell_line/"$selected_cell_line".fa -b $WD/"$condition_per_cell_line"/alignments_to_assembly/"$selected_cell_line"_nanocompore.bam --g $TRANSCRIPTOME_ASSEMBLY --samples --print-read-names --scale-events --iop $PARALLEL_JOBS > $WD_MODS/$selected_cell_line/nanocompore/alignments_to_assembly/eventalign/$condition_per_cell_line/events.tsv
        $NANOPOLISHCOMP NanopolishComp Eventalign_collapse -t $THREADS -i $WD_MODS/$selected_cell_line/nanocompore/alignments_to_assembly/eventalign/$condition_per_cell_line/events.tsv -o $WD_MODS/$selected_cell_line/nanocompore/alignments_to_assembly/eventalign/$condition_per_cell_line/collapse/

        # compare eventalign file for each cell line to IVT eventalign file, transcript per transcript
        for filename in $IVT_eventalign/eventalign*; do
                base=${filename##*/eventalign_}
                $NANOCOMPORE nanocompore sampcomp --file_list1 $WD_MODS/$selected_cell_line/nanocompore/alignments_to_assembly/eventalign/$condition_per_cell_line/collapse/out_eventalign_collapse.tsv  --file_list2 "$filename"/out_eventalign_collapse.tsv  --label1 $condition_per_cell_line --label2 IVT --fasta $TRANSCRIPTOME_ASSEMBLY --outpath $WD_MODS/$selected_cell_line/nanocompore/alignments_to_assembly/sampcomp/"$condition_per_cell_line"_vs_IVT/"$base" --overwrite --downsample_high_coverage 5000 --allow_warnings --min_coverage 30 --logit --nthreads $THREADS --bed $NANOCOMP_BED

        done
done

NANOCOMP_BED="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus_name_commas.bed"

