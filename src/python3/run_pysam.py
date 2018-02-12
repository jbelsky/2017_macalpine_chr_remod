#!/usr/bin/python3

import pandas as pd
import pysam
import os

# Set the files
f_bam = "/data/illumina_pipeline/aligned_experiments/DM278/DM278.bam"
f_acs = "/home/jab112/2017_macalpine_chr_remod/GenomeFeatureFiles/repliOriACS_sacCer2_OriDB_JABcurated_798sites.csv"
f_out = "/home/jab112/2017_macalpine_chr_remod/data/WT/DM278/DM278_WT_BrdU_OriDB_ACS798ori_1000bpWin.tsv"

# Set the window width
readWin = 1000

# Open the files
bam_pys = pysam.AlignmentFile(f_bam, mode = "rb")
acs_df = pd.read_csv(f_acs, sep = ",", header = 0, index_col = 0)

# Get the overall read depth
print("Obtaining total read counts across BAM...")
totalReads = bam_pys.count()
print("\t%s has %d total read counts" % (os.path.basename(f_bam), totalReads))

# Initialize the output df
output_df = acs_df.iloc[:,0:3]
output_df["rpkm"] = 0

acs_ct = 0

for f in acs_df.iterrows():

	acs_ct += 1
	if acs_ct % 25 == 0:
		print("Processing OriDB ACS:\t\t%d" % acs_ct, end = "\r")

	# Get the coordinates
	acs_name = f[0]
	chrom = str(f[1]["chr"])
	pos = f[1]["pos"]
	strand = f[1]["strand"]

	# Get the start and end coordinates
	l_pos = pos if pos - readWin > 0 else 0
	r_pos = pos if pos + readWin > bam_pys.get_reference_length(chrom) else bam_pys.get_reference_length(chrom)

	# Initialize the read_count
	read_ct = 0

	for alignRead in bam_pys.fetch(chrom, l_pos, r_pos):
		read_ct += 1

	# Obtain the RPKM
	output_df.loc[acs_name, "rpkm"] = read_ct / ((totalReads / 1E6) * ((r_pos - l_pos + 1) / 1000))

print("\nWriting to %s..." % f_out)
output_df.to_csv(f_out, sep = "\t", header = True, index = True, index_label = True, float_format = "%.4f")
print("\tComplete!")
