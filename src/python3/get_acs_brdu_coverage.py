#!/usr/bin/python3

import pandas as pd
import pysam
import os
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("design_file", help = "the chromatin remodeler design file")
parser.add_argument("acs_feature_file", help = "the ACS feature file")
parser.add_argument("--win", help = "+/- bp around feature position")
parser.add_argument("--out_dir", help = "output_directory")
args = parser.parse_args()

design_df = pd.read_table(args.design_file, sep = "\t", header = 0, index_col = "ID")

# Set the window width
readWin = int(args.win)

# Iterate through the double mutants
for dm_id in design_df[design_df["Mutant"] == "double"].index:

	print(dm_id + "\t" + design_df.loc[dm_id, "Genotype"])

	readWin = int(args.win)
	f_out = "/".join([args.out_dir, design_df.loc[dm_id,"Genotype"].replace(" ", "_"), dm_id, dm_id + "_brdu_coverage_2500bp_win.txt"])

	f_bam = design_df.loc[dm_id, "BamFile"]

	# Read in the chr positions
	pos_df = pd.read_table(args.acs_feature_file, sep = ",", header = 0, index_col = "name")
	pos_df["Counts"] = 0

	# Open the bam file
	bam_pys = pysam.AlignmentFile(f_bam, mode = "rb")

	for r_idx in pos_df.index:

		# Get the coordinates
		chrom = str(pos_df.loc[r_idx, "chr"])
		pos = pos_df.loc[r_idx, "pos"] 

		# Get the start and end coordinates of the acs
		l_pos = max(pos - readWin, 0)
		r_pos = min(pos + readWin, bam_pys.get_reference_length(chrom))

		# Initialize the read_count
		read_ct = 0

		for alignRead in bam_pys.fetch(chrom, max(l_pos - 1000, 0), min(r_pos + 1000, bam_pys.get_reference_length(chrom))):

			# Shift by 75 bp
			if alignRead.is_reverse:
				alignRead_pos = alignRead.pos - 75
			else:
				alignRead_pos = alignRead.pos + 75

			if alignRead_pos >= l_pos and alignRead_pos <= r_pos:
				read_ct += 1

		# Enter into the data frame
		pos_df.loc[r_idx, "Counts"] = read_ct

	# Write the output
	pos_df.to_csv(f_out, sep = "\t", header = True, index = True, columns = ["chr", "pos", "strand", "Counts"])
