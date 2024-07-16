#!/usr/bin/env python3

import os
import tskit 
import collections
import numpy as np
import allel

select_arr = [0.1, 0.1]
sigma_arr = [0.015, 0.5]


DIR = "/home2/mnc42/"
SUBDIR = "spatial_WF_grid_sweep_500k"
DOMINANCE = 0.5
SUFFIX = "_sampling_rd"
OUTPREFIX = DIR + SUBDIR + "/out"

window = 50000

samples = 100
genome_length = 1e7
total_files = 1000

tothap = 100
hc_arr = []
for select_idx, SELECTION in enumerate(select_arr):
	DISPERSAL_DISTANCE = sigma_arr[select_idx]
	outfile = OUTPREFIX + "_"  + str(SELECTION) + "_" + str(DISPERSAL_DISTANCE) + SUFFIX + ".hfs"
	output = open(outfile, "w")
	
	file_prefix =  DIR + SUBDIR + "/" + str(DISPERSAL_DISTANCE) + "_" + str(SELECTION) + "_" + str(DOMINANCE) + "_sampling"   
	files = os.listdir(file_prefix)
	files = [file for file in files if file.endswith(SUFFIX + ".trees")]
	
	print("total files = ", len(files))
	print(files[0])
	
	print(window)
	start = genome_length/2 - window/2
	stop = genome_length/2 + window/2

	allhc = []
	for file in files[:total_files]:
		nts = tskit.load(file_prefix + "/" + file)
		samps = np.random.choice(nts.samples(), samples, replace = False)
		subts = nts.keep_intervals([[start,stop]])
		h = allel.HaplotypeArray(subts.genotype_matrix()[:,:samples]) # get the haplotype matrix
		hc = h.distinct_counts()
		hc = list(hc[:tothap]) # only use first 10 haplotypes
		hc.extend([0 for i in range(tothap - len(hc))])# changed nan => 0 
		allhc.append(hc)
		output.write(str(DISPERSAL_DISTANCE) + "\t" + str(SELECTION) + "\t" + str(window) + "\t" + '\t'.join("{:.1f}".format(item) for item in hc) + "\n")
	output.close()
    

