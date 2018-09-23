#!/usr/bin/python3

import pandas as pd
import pysam
import os
import glob




# Set the window width
readWin = 3000

# Read in the chr positions
pos_df = pd.read_table("../../GenomeFeatureFiles/simulate_chr_pos.txt", sep = "\t", header = 0, index_col = None)

# Set the files
f_bam = "/data/illumina_pipeline/aligned_experiments/DM409/DM409.bam"

# Open the bam file
bam_pys = pysam.AlignmentFile(f_bam, mode = "rb")

# Initialize output list
output_read_ct = []

for r_idx, chrPos in pos_df.iterrows():

	# Get the coordinates
	chrom = str(chrPos["Chromosome"])
	pos = chrPos["Position"]

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
	output_read_ct.append(str(read_ct))


# Write the output
with open("out.txt", mode = "w") as fO:
	fO.write("\n".join(output_read_ct) + "\n")
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
