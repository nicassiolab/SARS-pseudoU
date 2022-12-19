#!/bin/bash
set -e -o pipefail
export LC_ALL=C

###      script to create .tar.gz files of the fast5 sequeunces for the in-house generated datasets

DRS_CaCo2_3="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210414_1407_X1_FAL82866_fb16777a/S33249_CaCo2_C34"
DRS_CaCo2_4="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210825_1506_X1_FAL77093_8c3d8e3b/S35755_CaCo2_C37"
WD="/hpcnfs/scratch/temporary/camilla_TL/SARS_CoV_2_ENA_data"

###	FAST5
#mkdir -p $WD/DRS_CaCo2_3
#cp -r $DRS_CaCo2_3/fast5_pass $WD/DRS_CaCo2_3
#cp -r $DRS_CaCo2_3/fast5_fail $WD/DRS_CaCo2_3
#tar -cvzf $WD/DRS_CaCo2_3.tar.gz $WD/DRS_CaCo2_3

#mkdir -p $WD/DRS_CaCo2_4
#cp -r $DRS_CaCo2_4/fast5_pass $WD/DRS_CaCo2_4
#cp -r $DRS_CaCo2_4/fast5_fail $WD/DRS_CaCo2_4
#tar -cvzf $WD/DRS_CaCo2_4.tar.gz $WD/DRS_CaCo2_4

###	FASTQ
cat $DRS_CaCo2_3/fastq_pass/*.fastq > $WD/DRS_CaCo2_3_pass.fastq
gzip $WD/DRS_CaCo2_3_pass.fastq

cat $DRS_CaCo2_4/fastq_pass/*.fastq > $WD/DRS_CaCo2_4_pass.fastq
gzip $WD/DRS_CaCo2_4_pass.fastq







