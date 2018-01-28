#!/bin/bash

# Set the DM_id
dm_id=( 593 594 596 597)

# Set the invariant parameters
feature_file_name="/data/genome_feature_files/eaton_acs_filter_rDNA_telomere_230_sites_ucsc_chr_names.csv"
jar_file_name="/opt/jar_files/chip-seq-analysis.jar"
program="ChIPSeqFeatureDensitySignal"
win=500
chip_shift=75
bw=20
span=1

for d in ${dm_id[@]}
do

	# Set the variables
	bam_file_name="/data/illumina_pipeline/aligned_experiments/DM${d}/DM${d}.bam"
	output_file_name="/data2/jab112/2015_yeast_chr_remod/mcm_chip_seq/mcm_chip_seq_DM${d}_around_eaton_acs_230_sites.csv"

	# Execute the script
	java -jar $jar_file_name $program $bam_file_name $feature_file_name $output_file_name $win $chip_shift $bw $span

done
