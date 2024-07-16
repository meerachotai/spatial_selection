#!/usr/bin/env python3

# PREFIX can be taken from .trees filename, sampled .trees file will be named similarly
# /home2/mnc42/sampling_strategies_inds_dest.py ${PREFIX} ${MU} ${SAMPLE_SIZE} "sampling" ${SEED} --tree -mid
# /home2/mnc42.py /home2/mnc42/ spatial_WF_grid_sweep_from_neutral /home2/mnc42/spatial_WF_grid_sweep_from_neutral/sfs_window -rd
# /home2/mnc42/SGV_mutation_age.py 0.1 /home2/mnc42/ SGV_pop_1 /home2/mnc42/SGV_pop_1/out

import argparse
import sys
import os
import numpy as np
import tskit
import msprime
import pyslim
import glob
from tqdm import tqdm
import random

parser = argparse.ArgumentParser()
parser.add_argument("selection", help="selection coef / 0 if neutral", type = float)
parser.add_argument("dir", help="main directory", type = str)
parser.add_argument("subdir", help="subdirectory", type = str)
parser.add_argument("out", help="outfile (including path)", type = str, default = "out")

args = parser.parse_args()
SELECTION = args.selection
DIR=args.dir
SUBDIR=args.subdir
OUT=args.out


# DIR = "/home2/mnc42/"
# SUBDIR = "spatial_WF_grid_sweep_from_neutral"

# 0.05 0.01 0.005 0.001 0.0005 0.0001
# disp_arr = [0.05, 0.01, 0.005, 0.001, 0.0001]
disp_arr = [0.03, 0.04]
DOMINANCE = 0.5

SAMPLE_SIZE = 100
LENGTH = 1e7

N = 1000 # number of replicates
SUFFIX = "sampling_rd"

for FREQUENCY in disp_arr:
	OUTFILE = OUT + "_" + str(FREQUENCY) + "_" + str(SELECTION) + "_" + str(DOMINANCE) + "_" + SUFFIX
	output = open(OUTFILE + ".age", "w")

	file_prefix = DIR + SUBDIR + "/"+ str(SELECTION) + "/" + str(SELECTION) + "_" + str(DOMINANCE) + "_" + str(FREQUENCY) + "*"+ str(SAMPLE_SIZE) + "_" + SUFFIX + "_*.trees"

	print(file_prefix)
	files = glob.glob(file_prefix)
	print(len(files))
# 	files = np.random.choice(files, N)
	random.shuffle(files)
	time = []
	slim_time = []
	
	for file in tqdm(files):
		ts = tskit.load(file)
		pos = file.split("_")[-1].split(".")[0]
		if(pos != "neutral" and pos != ""):
			time.append(ts.site(position = int(pos)).mutations[0].time)
			slim_time.append(ts.site(position = int(pos)).mutations[0].metadata['mutation_list'][0]['slim_time'])
	print(len(time))
	res = [str(int(x)) + " " + str(int(-y)) for x, y in zip(time, slim_time)]
	output.write('\n'.join(item for item in res) + "\n")
	
	output.close()