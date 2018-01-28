#!/bin/bash

# Set the variables
bam_file_name="/data/illumina_pipeline/projects/rachel_dros_mnase/asynch/dros_mnase_asynch_DM464_-m1.bam"
feature_file_name="/data/illumina_pipeline/scripts/feature_files/dros/orc_chip_seq/kc_orc_chip_seq_peak.csv"
output_file_name="/data2/jab112/2015_rachel_dros_mnase/DM464_-m1_0_120_frag_bw-10_around_kc_orc_chip_seq_peak.csv"

# Set the parameters
win=500
fragL=0
fragH=120
bw=10

# Set the jar file path
jar_file_path="/opt/jar_files/java-mnase-seq-analysis.jar"

# Execute the script
java -jar $jar_file_path $bam_file_name $feature_file_name $output_file_name $win $fragL $fragH $bw
