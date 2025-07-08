#!/usr/bin/env python3

# /home2/mnc42/sampling_strategies_SGV_present_inds.py $PREFIX $MU $SAMPLE_SIZE "sampling" $SEED --tree -rd
# /home2/mnc42/haplotype_stats_SGV.py /home2/mnc42/ SGV_pop_1 /home2/mnc42/SGV_pop_1/out -rd

import argparse
import sys
import os
import numpy as np
import tskit
import glob
from tqdm import tqdm
import allel
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("dir", help="main directory", type = str)
parser.add_argument("subdir", help="subdirectory", type = str)
parser.add_argument("selection", help="selection coef", type = float)
parser.add_argument("out", help="outfile (including path)", type = str, default = "out")
parser.add_argument("nfiles", help="number of files", type = int)
parser.add_argument("nsites", help="number of sites within window", type = int, default = 1000)

args = parser.parse_args()
DIR=args.dir
SUBDIR=args.subdir
OUT=args.out
N = args.nfiles # number of replicates
SELECTION = args.selection
NSITES = args.nsites

freq_arr = ["0.00005", "0.05"] # starting frequency
SELECTION = 0.1

SAMPLE_SIZE = 100 # number of samples

DOMINANCE = 0.5
SAMPLE_SIZE = 100
genome_length = 1e7

SUFFIX = "sampling_rd"	

# OUTFILE = OUT + "_" + str(N) + "_" + str(SAMPLE_SIZE) + "_" + str(NSITES) + "_" + SUFFIX
OUTFILE = OUT + "_" + str(SELECTION) + "_" + str(N) + "_" + str(SAMPLE_SIZE) + "_" + str(NSITES) + "_" + SUFFIX

output_hstats = open(OUTFILE + "_SGV_pop.garudh", "w")


for freq_idx,FREQUENCY in enumerate(freq_arr):

	file_prefix = DIR + SUBDIR + "/"+ str(SELECTION) + "/" + str(SELECTION) + "_" + str(DOMINANCE) + "_" + str(FREQUENCY) + "*"+ str(SAMPLE_SIZE) + "_" + SUFFIX + "_*.trees"
	print(file_prefix)

	files = glob.glob(file_prefix)  # use files generated from sampling_strategies_inds.py
	print(len(files)); 
	files = np.random.choice(files, N,replace = False)


	for file in tqdm(files):
		nts = tskit.load(file)
		selected_site = int(file.split("_")[-1].split(".")[0])
		
		if(file == files[0]):
			print(selected_site)
			
		h = allel.HaplotypeArray(nts.genotype_matrix())
		pos = nts.sites_position

		idx = np.argwhere(pos == selected_site)[0][0]
		START = idx - (NSITES//2)
		STOP = idx + (NSITES//2)
		
		hstat = allel.moving_garud_h(h,size = NSITES, start = START, stop = STOP)
		
		output_hstats.write(str(FREQUENCY) + "\t" + str(SELECTION) + "\t" + str(NSITES) + "\t" + '\t'.join("{:.10f}".format(item[0]) for item in hstat) + "\n")


output_hstats.close()


