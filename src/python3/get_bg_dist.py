#!/usr/bin/python3

import pandas as pd
import numpy as np
import pysam
import os
import glob
import random

dmID = "515"
f_bam = "/data/illumina_pipeline/aligned_experiments/DM" + dmID + "/dm" + dmID + ".bam"

bam_pys = pysam.AlignmentFile(f_bam, mode = "rb")

chr_ref = np.array(bam_pys.header.references)
chr_lengths = np.array(bam_pys.header.lengths)
chr2len_dict = dict(zip(chr_ref, chr_lengths))

# Get the probability distribution of the chromosomes
chr_prob = chr_lengths / np.sum(chr_lengths)

sim = np.random.choice(chr_ref, size = 1000, replace = True, p = chr_prob)

# Get the counts of each chromosome
ref, counts = np.unique(sim, return_counts = True)
ref_counts = dict(zip(ref, counts))

sim_pos = {}

for r,c in ref_counts.items():
	print(r + "\t" + str(c) + "\t" + str(chr2len_dict[r]))

	pos = random.sample(range(3001, chr2len_dict[r] - 3000), c)
	pos.sort()
	sim_pos[int(r)] = pos

# Get the chromosomes sorted
chr_ref_int = list(sim_pos.keys())
chr_ref_int.sort()

# Write the output
with open("../../GenomeFeatureFiles/simulate_chr_pos.txt", mode = "w") as fO:
	for c in chr_ref_int:
		for p in sim_pos[c]:
			fO.write(str(c) + "\t" + str(p) + "\n")

fO.close()

'''

dm_id = [510, 511, 512, 513, 514, 516, 517, 518]
chr_remod = ["mrc1_isw1_isw2_chd1", "mrc1_isw1_isw2_asf1", "mrc1_isw1_chd1_asf1", "mrc1_isw2_swr1_asf1",
			  "mrc1_isw2_swr1_isw1", "mrc1_chd1_swr1_isw1", "mrc1_swr1_asf1_isw1", "mrc1_isw2_chd1_swr1"
			]

for idx in range(0, len(dm_id)):

	dmID = str(dm_id[idx])
	remod = chr_remod[idx]

	# Set the files
	f_bam = "/data/illumina_pipeline/aligned_experiments/DM" + dmID + "/dm" + dmID + ".bam"
	f_acs = "/home/jab112/2017_macalpine_chr_remod/GenomeFeatureFiles/repliOriACS_sacCer2_OriDB_JABcurated_798sites.csv"
	f_out = "/home/jab112/2017_macalpine_chr_remod/data/" + remod + "/DM" + dmID + "/DM" + dmID + "_" + remod + "_BrdU_OriDB_ACS798ori_3000bpWin.tsv"

	# Set the window width
	readWin = 3000

	#exp_id = "DM273"
	#mutant = "single"
	#mrc1 = "wt"
	#genotype = "WT"
	#replicate = 1

	d_bam = "/data/illumina_pipeline/aligned_experiments"
	p = re.complie("DM273", re.IGNORECASE)

	fInDir = os.listdir(d_bam + "/" + exp_id)
	f_bam = ""
	for f in fInDir:
		if p.match(f):
			f_bam = f

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

		# Get the start and end coordinates of the acs
		l_pos = max(pos - readWin, 0)
		r_pos = min(pos + readWin, bam_pys.get_reference_length(chrom))

		# Get the buffer search region
		l_pos_search = l_pos - 1000
		r_pos_search = r_pos + 1000

		# Initialize the read_count
		read_ct = 0

		for alignRead in bam_pys.fetch(chrom, max(l_pos_search, 0), min(r_pos_search, bam_pys.get_reference_length(chrom))):
			# Shift by 75 bp
			if alignRead.is_reverse:
				alignRead_pos = alignRead.pos - 75
			else:
				alignRead_pos = alignRead.pos + 75

			if alignRead_pos >= l_pos and alignRead_pos <= r_pos:
				read_ct += 1

		# Obtain the RPKM
		output_df.loc[acs_name, "rpkm"] = read_ct / ((totalReads / 1E6) * ((r_pos - l_pos + 1) / 1000))

	print("\nWriting to %s..." % f_out)
	output_df.to_csv(f_out, sep = "\t", header = True, index = True, index_label = "Samples", float_format = "%.4f")
	print("\tComplete!")


'''
