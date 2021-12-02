#!/bin/bash
set -e -o pipefail
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
TRANSCRIPTOME_ASSEMBLY="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
TOOLS="/hpcnfs/scratch/FN/TL/cugolini/tools"
MISEQ="/hpcnfs/techunits/genomics/PublicData/Leonardi/tleonardi/FASTQ/211021_M03114_0038_000000000-DDBBR/Sample_S37211_ORF_9D_PCR_product"

###	trim the 5' and 3' adapter
#mkdir -p $WD/miseq/FASTQ
#$TOOLS/reaper -i $MISEQ/S37211_ORF_9D_PCR_product_S1_L001_R1_001.fastq.gz -geom no-bc -3pa GATCGGAAGAGCACACGTC -basename $WD/miseq/FASTQ/R1
#$TOOLS/reaper -i $MISEQ/S37211_ORF_9D_PCR_product_S1_L001_R2_001.fastq.gz -geom no-bc -3pa CGGTGGTCGCCGTATCATT -basename $WD/miseq/FASTQ/R2
#sed -n '1~4s/^@/>/p;2~4p' $WD/miseq/FASTQ/R1.lane.clean > $WD/miseq/FASTQ/R1.fasta		### from fastq to fasta
#sed -n '1~4s/^@/>/p;2~4p' $WD/miseq/FASTQ/R2.lane.clean > $WD/miseq/FASTQ/R2.fasta
#awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' $WD/miseq/FASTQ/R1.fasta > $WD/miseq/FASTQ/R1_filt.fasta		### rm 0-length sequences
#awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' $WD/miseq/FASTQ/R2.fasta > $WD/miseq/FASTQ/R2_filt.fasta
#rm $WD/miseq/FASTQ/R1.fasta
#rm $WD/miseq/FASTQ/R2.fasta

#mkdir -p $WD/miseq/map_to_genome/ALIGNMENTS_BACKUP
#mkdir -p  $DATA/SARS_CoV_2_genome/indices
#cp $GENOME_FA $DATA/SARS_CoV_2_genome/
#singularity pull docker://tleonardi/star:v2.9.7a
#mv star_v2.9.7a.sif /hpcnfs/scratch/TSSM/cugolini/cov/img/
#star="/hpcnfs/scratch/TSSM/cugolini/cov/img/star_v2.9.7a.sif"
#singularity exec -B /hpcnfs/scratch/ $star STAR --runMode genomeGenerate --genomeDir $DATA/SARS_CoV_2_genome/indices --genomeFastaFiles $DATA/SARS_CoV_2_genome/edited.fa --genomeSAindexNbases 7
#singularity exec -B /hpcnfs/scratch/ $star STAR --runThreadN 8 --genomeDir $DATA/SARS_CoV_2_genome/indices --readFilesIn $WD/miseq/FASTQ/R1_filt.fasta $WD/miseq/FASTQ/R2_filt.fasta --outFilterIntronStrands None --outSJfilterOverhangMin 10 12 12 12 --outFileNamePrefix $WD/miseq/map_to_genome/
#source activate $ENVS/nanopore_pipeline_env/
#samtools view -h -Sb $WD/miseq/map_to_genome/Aligned.out.sam > $WD/miseq/map_to_genome/ALIGNMENTS_BACKUP/Aligned.out.bam
#samtools view -F 2316 -Sb $WD/miseq/map_to_genome/Aligned.out.sam > $WD/miseq/map_to_genome/Aligned.out.bam
#bedtools bamtobed -bed12 -i $WD/miseq/map_to_genome/Aligned.out.bam > $WD/miseq/map_to_genome/Aligned.out.bed
#awk 'BEGIN{ORS=" "}{if($10<4){split($12,a,","); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11;for(i=1;i<=length(a);i++) print a[i]+$2; print "\n"}}' $WD/miseq/map_to_genome/Aligned.out.bed| sort -k10,10n -k12,12n -k13,13n -k14,14n | awk -v OFS="\t" '$1=$1' | awk 'BEGIN{ORS=" "}{print $1,$2,$3,$4,$5,"+",$7,$8,$9,$10,$11; for(i=12;i<=NF;i++) printf "%s," ,($i-$2); print "\n"}' |awk -v OFS="\t" '$1=$1'| awk 'BEGIN{OFS="\t"}{$12=substr($12,1,length($12)-1); print $0}' | awk '$10>1' > $WD/miseq/map_to_genome/junction_sorted_2_3_ex.bed	###awk script to sort the bed file according to junctions
#bedtools genomecov -i $WD/miseq/map_to_genome/twoex_track_canonical.bed -g genomefile.genome -bg -trackline -strand + > $WD/miseq/map_to_genome/twoex_track_canonical.bedgraph


###	create .tar.gz file to submit
mkdir -p /hpcnfs/scratch/temporary/camilla_TL/SARS_CoV_2_ENA_data/MiSeq/
#cp $MISEQ/S37211_ORF_9D_PCR_product_S1_L001_R1_001.fastq.gz  /hpcnfs/scratch/temporary/camilla_TL/SARS_CoV_2_ENA_data/MiSeq/
#cp $MISEQ/S37211_ORF_9D_PCR_product_S1_L001_R2_001.fastq.gz  /hpcnfs/scratch/temporary/camilla_TL/SARS_CoV_2_ENA_data/MiSeq/
#gunzip  /hpcnfs/scratch/temporary/camilla_TL/SARS_CoV_2_ENA_data/MiSeq/S37211_ORF_9D_PCR_product_S1_L001_R1_001.fastq.gz 
#gunzip /hpcnfs/scratch/temporary/camilla_TL/SARS_CoV_2_ENA_data/MiSeq/S37211_ORF_9D_PCR_product_S1_L001_R2_001.fastq.gz
tar -czvf /hpcnfs/scratch/temporary/camilla_TL/SARS_CoV_2_ENA_data/MiSeq.tar.gz /hpcnfs/scratch/temporary/camilla_TL/SARS_CoV_2_ENA_data/MiSeq

