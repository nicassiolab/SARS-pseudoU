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
singpore="/hpcnfs/scratch/TSSM/cugolini/tools/porechop/porechop.simg"

###	SAMPLES

#MATTHEWS="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/eventalign/eventalign_matthews_caco/out_eventalign_collapse.tsv"
#SRAFF="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/eventalign/eventalign_sraf_caco2/out_eventalign_collapse.tsv"
#WT_C37="/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/eventalign/WT/out_eventalign_collapse.tsv"
#WT_C34="/hpcnfs/scratch/FN/TL/cugolini/cov/analysis/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/WT/out_eventalign_collapse.tsv"

#mkdir -p $WD/NANOCOMPORE/sampcomp/caco2

#for filename in /hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/eventalign/IVT/eventalign*; do
#       base=${filename##*/eventalign_}
#       singularity exec -B /hpcnfs/scratch/ $nanocomp104 nanocompore sampcomp --file_list1 $WT_C34,$WT_C37,$SRAFF,$MATTHEWS --file_list2 "$filename"/out_eventalign_collapse.tsv  --label1 WT --label2 IVT --fasta $NANOCOMP_FA --outpath $WD/NANOCOMPORE/sampcomp/caco2/"$base" --overwrite --downsample_high_coverage 5000 --allow_warnings --min_coverage 30 --logit --nthreads 10 --bed $NANOCOMP_BED
#done

#while read filename; do
#        base=${filename##*/eventalign_}
#        singularity exec -B /hpcnfs/scratch/ $nanocomp104 nanocompore sampcomp --file_list1 $MATTHEWS,$SRAFF,$WT_C37,$WT_C34 --file_list2 "$filename"/out_eventalign_collapse.tsv  --label1 WT --label2 IVT --fasta $NANOCOMP_FA --outpath $WD/NANOCOMPORE/sampcomp/caco2/"$base" --overwrite --downsample_high_coverage 5000 --allow_warnings --min_coverage 30 --logit --nthreads 11 --bed $NANOCOMP_BED
#done</hpcnfs/scratch/FN/TL/cugolini/cov/scripts/caco2_tx.txt


#############################################	CAT ALL SAMPLES TOGETHER


# map to the assembly(EXTRACTION)
mkdir -p $WD/caco2/map_to_recap_assembly

#cat /hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20200703_1147_X2_FAL55849_d6d1028f/S29022_Caco_2_INF_MOI_0_1/fastq_pass/*.fastq /hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210414_1407_X1_FAL82866_fb16777a/S33249_CaCo2_C34/fastq_pass/*.fastq /hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210825_1506_X1_FAL77093_8c3d8e3b/S35755_CaCo2_C37/fastq_pass/*.fastq | sed -n '1~4s/^@/>/p;2~4p' > $WD/caco2/map_to_recap_assembly/temp.fa
#cat /hpcnfs/scratch/TSSM/tleonardi/SARS-CoV-2-datasets/Matthews_liverpool_caco/fastq/ont_research.fa $WD/caco2/map_to_recap_assembly/temp.fa > $WD/caco2/map_to_recap_assembly/caco2.fa
#rm $WD/caco2/map_to_recap_assembly/temp.fa
#singularity exec -B /hpcnfs/scratch/ $singpore minimap2 -t 20 -ax map-ont -p 0 -N 10 $TRANSCRIPTOME_ASSEMBLY $WD/caco2/map_to_recap_assembly/caco2.fa > $WD/caco2/map_to_recap_assembly/caco2.sam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -Sb $WD/caco2/map_to_recap_assembly/caco2.sam > $WD/caco2/map_to_recap_assembly/caco2.bam
#rm $WD/caco2/map_to_recap_assembly/caco2.sam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -F 2324 -Sb $WD/caco2/map_to_recap_assembly/caco2.bam | singularity exec -B /hpcnfs/scratch/ $singpore samtools sort > $WD/caco2/map_to_recap_assembly/caco2_filtered.bam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools index $WD/caco2/map_to_recap_assembly/caco2_filtered.bam

### F5C (eventalign)

#mkdir -p $WD/caco2/temporary_fast5/
mkdir -p  $WD/caco2/NANOCOMPORE/eventalign/collapse/
mkdir -p  $WD/caco2/NANOCOMPORE/sampcomp/
FAST5_TOTAL="$WD/caco2/temporary_fast5/"
PUBLIC_DATA="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi"

#cp /hpcnfs/scratch/TSSM/tleonardi/SARS-CoV-2-datasets/Matthews_liverpool_caco/fast5/*.fast5 $FAST5_TOTAL
#cp $PUBLIC_DATA/FAST5/20210825_1506_X1_FAL77093_8c3d8e3b/S35755_CaCo2_C37/fast5_pass/*.fast5 $FAST5_TOTAL
#cp $PUBLIC_DATA/FAST5/20200703_1147_X2_FAL55849_d6d1028f/S29022_Caco_2_INF_MOI_0_1/fast5_pass/*.fast5 $FAST5_TOTAL
#cp $PUBLIC_DATA/FAST5/20210414_1407_X1_FAL82866_fb16777a/S33249_CaCo2_C34/fast5_pass/*.fast5 $FAST5_TOTAL

#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux index -t 20 -d $FAST5_TOTAL $WD/caco2/map_to_recap_assembly/caco2.fa
#/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux eventalign --rna --min-mapq 0 -t 10 -r $WD/caco2/map_to_recap_assembly/caco2.fa -b $WD/caco2/map_to_recap_assembly/caco2_filtered.bam --g $TRANSCRIPTOME_ASSEMBLY --samples --print-read-names --scale-events --iop 5 > $WD/caco2/NANOCOMPORE/eventalign/events.tsv
singularity exec -B /hpcnfs/scratch/ /hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_8d6a70c.img  NanopolishComp Eventalign_collapse -t 10 -i $WD/caco2/NANOCOMPORE/eventalign/events.tsv -o $WD/caco2/NANOCOMPORE/eventalign/collapse/

for filename in /hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/eventalign/IVT/eventalign*; do
       base=${filename##*/eventalign_}
       singularity exec -B /hpcnfs/scratch/ $nanocomp104 nanocompore sampcomp --file_list1 $WD/caco2/NANOCOMPORE/eventalign/collapse/out_eventalign_collapse.tsv  --file_list2 "$filename"/out_eventalign_collapse.tsv  --label1 WT --label2 IVT --fasta $NANOCOMP_FA --outpath  $WD/caco2/NANOCOMPORE/sampcomp/"$base" --overwrite --downsample_high_coverage 5000 --allow_warnings --min_coverage 30 --logit --nthreads 10 --bed $NANOCOMP_BED
done




