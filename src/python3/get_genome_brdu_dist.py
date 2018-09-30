#!/usr/bin/python3

import pandas as pd
import pysam
import os
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("design_file", help = "the chromatin remodeler design file")
parser.add_argument("sim_genome_pos_file", help = "simulated genome positions")
parser.add_argument("--win", help = "+/- bp around feature position")
parser.add_argument("--out_dir", help = "output_directory")
args = parser.parse_args()

design_df = pd.read_table(args.design_file, sep = "\t", header = 0, index_col = "ID")

dm_id = "DM409"

# Set the window width
readWin = int(args.win)
f_out = "/".join([args.out_dir, design_df.loc[dm_id,"Genotype"], dm_id, dm_id + "_sim.txt"])

f_bam = design_df.loc[dm_id, "BamFile"]


# Read in the chr positions
pos_df = pd.read_table(args.sim_genome_pos_file, sep = "\t", header = 0, index_col = None)
pos_df["Counts"] = 0

# Open the bam file
bam_pys = pysam.AlignmentFile(f_bam, mode = "rb")


for r_idx in pos_df.index:

	# Get the coordinates
	chrom = str(pos_df.loc[r_idx, "Chromosome"])
	pos = pos_df.loc[r_idx, "Position"]

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

	# Enter into output list
	pos_df.loc[r_idx, "Counts"] = read_ct

	if r_idx % 500 == 0:
		print("Processed %d simulated positions..." % r_idx, end = "\r")

# Write the output
pos_df.to_csv(f_out, sep = "\t", header = True, index = False)

print("\n\n\tComplete!")
