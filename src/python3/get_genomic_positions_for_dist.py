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

sim = np.random.choice(chr_ref, size = 10000, replace = True, p = chr_prob)

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
	fO.write("Chromosome\tPosition\n")
	for c in chr_ref_int:
		for p in sim_pos[c]:
			fO.write(str(c) + "\t" + str(p) + "\n")

fO.close()
