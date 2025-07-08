#!/usr/bin/env python3

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
parser.add_argument("-m","--mid", help="midpoint chromosome sampling", action = "store_true")
parser.add_argument("-rd","--random_diploid", help="random individual sampling", action = "store_true")

args = parser.parse_args()
DIR=args.dir
SUBDIR=args.subdir
OUT=args.out
N = args.nfiles # number of replicates
SELECTION = args.selection
NSITES = args.nsites

# disp_arr = [0.015, 0.02, 0.04, 0.08, 0.1, 0.5]
disp_arr = [0.015, 0.5]

SAMPLE_SIZE = 100 # number of samples

DOMINANCE = 0.5
SAMPLE_SIZE = 100
genome_length = 1e7
selected_site = genome_length // 2


if(args.mid):
	SUFFIX = "sampling_mid"
if(args.random_diploid):
	SUFFIX = "sampling_rd"	# for neutral we would want sampling_*_rd

OUTFILE = OUT + "_" + str(SELECTION) + "_" + str(N) + "_" + str(SAMPLE_SIZE) + "_" + str(NSITES) + "_" + SUFFIX
output_hstats = open(OUTFILE + ".garudh", "w")

for disp_idx,DISPERSAL_DISTANCE in enumerate(disp_arr):
	file_prefix = DIR + SUBDIR + "/" + str(DISPERSAL_DISTANCE) + "_" + str(SELECTION) + "_" + str(DOMINANCE) + "_sampling" + "/*"+ str(SAMPLE_SIZE) + "_" + SUFFIX + ".trees"
	
	files = glob.glob(file_prefix)  # use files generated from sampling_strategies_inds.py
	files = np.random.choice(files, N, replace = False)

	print(len(files)); print(file_prefix)
	
	for file in tqdm(files):
		nts = tskit.load(file)
		
		h = allel.HaplotypeArray(nts.genotype_matrix())
		pos = nts.sites_position
		
		idx = np.argwhere(pos == selected_site)[0][0]
		START = idx - (NSITES//2)
		STOP = idx + (NSITES//2)
		
		hstat = allel.moving_garud_h(h,size = NSITES, start = START, stop = STOP)
		
		output_hstats.write(str(DISPERSAL_DISTANCE) + "\t" + str(SELECTION) + "\t" + str(NSITES) + "\t" + '\t'.join("{:.10f}".format(item[0]) for item in hstat) + "\n")
		

output_hstats.close()


