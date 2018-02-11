#!/bin/bash

# Set the invariant parameters
feature_file_name="/data/git/2017_macalpine_chr_remod/GenomeFeatureFiles/repliOriACS_sacCer2_OriDB_JABcurated_798sites.csv"
jar_file_name="/data/java-library/chip-seq-analysis.jar"
program="ChIPSeqFeatureDensitySignal"
win=5000
chip_shift=75
bw=50
span=10

for d in 409
do

	# Create the output directory
	mkdir -p -v "/data/git/2017_macalpine_chr_remod/data/mrc1/DM409"

	# Set the variables
	bam_file_name="/data/illumina_pipeline/aligned_experiments/DM${d}/DM${d}.bam"
	output_file_name="/data/git/2017_macalpine_chr_remod/data/mrc1/DM${d}/DM${d}_mrc1_BrdU_OriDB_ACS798ori.csv"

	# Execute the script
	java -jar $jar_file_name $program $bam_file_name $feature_file_name $output_file_name $win $chip_shift $bw $span

done
