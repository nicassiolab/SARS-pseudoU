#!/bin/bash
set -e -o pipefail
#export LC_ALL=C

source /hpcnfs/home/ieo5215/miniconda/etc/profile.d/conda.sh

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
DATA="$BASEDIR/data"
WD="$BASEDIR/analysis/per_cell_line"
REF_DATA="/hpcnfs/scratch/FN/camilla/nanopore/data"
ENVS="/hpcnfs/home/ieo5215/miniconda/envs"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
TRANSCRIPTOME_ASSEMBLY="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
NANOCOMP_FA="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/consensus_extracted.fa"
NANOCOMP_BED="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus_name_commas.bed"
nanocomp104="/hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_v1.0.4.sif"
#singularity pull tleonardi/nanocompore:v1.0.3
singpore="/hpcnfs/scratch/TSSM/cugolini/tools/porechop/porechop.simg"
###	SAMPLES

DAVIDSON="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/eventalign/eventalign_davidson/out_eventalign_collapse.tsv"
KIM="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/eventalign/eventalign_kim/out_eventalign_collapse.tsv"
TAIAROA="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/eventalign/eventalign_taiaroa/out_eventalign_collapse.tsv"
MATTHEWS="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/eventalign/eventalign_matthews_vero/out_eventalign_collapse.tsv"

mkdir -p $WD/NANOCOMPORE/sampcomp/vero

#while read filename; do
#	base=${filename##*/eventalign_}
#	singularity exec -B /hpcnfs/scratch/ $nanocomp104 nanocompore sampcomp --file_list1 $DAVIDSON,$KIM,$TAIAROA,$MATTHEWS --file_list2 "$filename"/out_eventalign_collapse.tsv  --label1 WT --label2 IVT --fasta $NANOCOMP_FA --outpath $WD/NANOCOMPORE/sampcomp/vero/"$base" --overwrite --downsample_high_coverage 5000 --allow_warnings --min_coverage 30 --logit --nthreads 12 --bed $NANOCOMP_BED
#done</hpcnfs/scratch/FN/TL/cugolini/cov/scripts/vero_tx.txt



###########   CAT ALL SAMPLES TOGETHER	##########


# map to the assembly(EXTRACTION)
mkdir -p $WD/vero/map_to_recap_assembly

DAVIDSON_FA="/hpcnfs/scratch/TSSM/tleonardi/SARS-CoV-2-datasets/Davidson_2020/fastq/ont_research.fa"
KIM_FA="/hpcnfs/scratch/TSSM/tleonardi/SARS-CoV-2-datasets/Kim_2020/fastq/Vero-Infected/ont_research.fa"
TAIAROA_FA="/hpcnfs/scratch/TSSM/tleonardi/SARS-CoV-2-datasets/Taiaroa_2020/fastq/CellCulture/ont_research.fa"
MATTHEWS_FA="/hpcnfs/scratch/TSSM/tleonardi/SARS-CoV-2-datasets/Matthews_liverpool_vero/fastq/ont_research.fa"


#cat $DAVIDSON_FA $KIM_FA $TAIAROA_FA $MATTHEWS_FA > $WD/vero/map_to_recap_assembly/vero.fa
#singularity exec -B /hpcnfs/scratch/ $singpore minimap2 -t 10 -ax map-ont -p 0 -N 10 $TRANSCRIPTOME_ASSEMBLY $WD/vero/map_to_recap_assembly/vero.fa > $WD/vero/map_to_recap_assembly/vero.sam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -Sb $WD/vero/map_to_recap_assembly/vero.sam > $WD/vero/map_to_recap_assembly/vero.bam
#rm $WD/vero/map_to_recap_assembly/vero.sam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -F 2324 -Sb $WD/vero/map_to_recap_assembly/vero.bam | singularity exec -B /hpcnfs/scratch/ $singpore samtools sort > $WD/vero/map_to_recap_assembly/vero_filtered.bam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools index $WD/vero/map_to_recap_assembly/vero_filtered.bam

### F5C (eventalign)

mkdir -p $WD/vero/temporary_fast5/
mkdir -p $WD/vero/NANOCOMPORE/eventalign/
FAST5_TOTAL="$WD/vero/temporary_fast5/"
PUBLIC_DATA="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi"
#cp /hpcnfs/scratch/TSSM/tleonardi/SARS-CoV-2-datasets/Matthews_liverpool_vero/fast5/*.fast5 $FAST5_TOTAL
#cp /hpcnfs/scratch/TSSM/tleonardi/SARS-CoV-2-datasets/Kim_2020/fast5_uncompressed/*.fast5 $FAST5_TOTAL
#cp /hpcnfs/scratch/TSSM/tleonardi/SARS-CoV-2-datasets/Taiaroa_2020/fast5/CellCulture/pass/*.fast5 $FAST5_TOTAL
#cp /hpcnfs/scratch/TSSM/tleonardi/SARS-CoV-2-datasets/Matthews_liverpool_vero/fast5/*.fast5 $FAST5_TOTAL
#chmod +xwr $FAST5_TOTAL

#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux index -t 10 -d $FAST5_TOTAL $WD/vero/map_to_recap_assembly/vero.fa
#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux eventalign --rna --min-mapq 0 -t 10 -r $WD/vero/map_to_recap_assembly/vero.fa -b $WD/vero/map_to_recap_assembly/vero_filtered.bam --g $TRANSCRIPTOME_ASSEMBLY --samples --print-read-names --scale-events --iop 5 > $WD/vero/NANOCOMPORE/eventalign/events.tsv
#singularity exec -B /hpcnfs/scratch/ /hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_8d6a70c.img  NanopolishComp Eventalign_collapse -t 10 -i $WD/vero/NANOCOMPORE/eventalign/events.tsv -o $WD/vero/NANOCOMPORE/eventalign/collapse/

for filename in /hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/eventalign/IVT/eventalign*; do
       base=${filename##*/eventalign_}
       singularity exec -B /hpcnfs/scratch/ $nanocomp104 nanocompore sampcomp --file_list1 $WD/vero/NANOCOMPORE/eventalign/collapse/out_eventalign_collapse.tsv  --file_list2 "$filename"/out_eventalign_collapse.tsv  --label1 WT --label2 IVT --fasta $NANOCOMP_FA --outpath  $WD/vero/NANOCOMPORE/sampcomp/"$base" --overwrite --downsample_high_coverage 5000 --allow_warnings --min_coverage 30 --logit --nthreads 10 --bed $NANOCOMP_BED
done



