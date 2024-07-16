#!/usr/bin/env python3

# PREFIX can be taken from .trees filename, sampled .trees file will be named similarly
# /home2/mnc42/sampling_strategies_inds.py ${PREFIX} ${MU} ${SAMPLE_SIZE} "sampling" ${SEED} --tree -rd
# /home2/mnc42/haplotype_stats_dispersal.py /home2/mnc42/ spatial_WF_grid_sweep_from_neutral /home2/mnc42/spatial_WF_grid_sweep_from_neutral/out -rd

import argparse
import sys
import os
import numpy as np
import tskit
import glob
from tqdm import tqdm
import collections
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("dir", help="main directory", type = str)
parser.add_argument("subdir", help="subdirectory", type = str)
parser.add_argument("selection", help="subdirectory", type = float)
parser.add_argument("out", help="outfile (including path)", type = str, default = "out")
parser.add_argument("nfiles", help="number of files", type = int)
parser.add_argument("-m","--mid", help="midpoint chromosome sampling", action = "store_true")
parser.add_argument("-rd","--random_diploid", help="random individual sampling", action = "store_true")

args = parser.parse_args()
DIR=args.dir
SUBDIR=args.subdir
OUT=args.out
N = args.nfiles # number of replicates
SELECTION = args.selection

disp_arr = [0.015, 0.02, 0.04, 0.08, 0.1, 0.5]
# select_arr = [0.01, 0.05, 0.1]
window_sizes = [10000, 30000,100000]

SAMPLE_SIZE = 100 # number of samples
# DIR = "/home2/mnc42/"
# SUBDIR = "spatial_WF_grid_sweep_from_neutral"

DOMINANCE = 0.5
SAMPLE_SIZE = 100
SEED = 5
genome_length = 1e7
n = 2

if(args.mid):
	SUFFIX = "sampling_mid"
if(args.random_diploid):
	SUFFIX = "sampling_rd"	# for neutral we would want sampling_*_rd

OUTFILE = OUT + "_" + str(SELECTION) + "_" + str(N) + "_" + str(SAMPLE_SIZE) + "_" + SUFFIX
output_gs = open(OUTFILE + ".het", "w")
output_nh = open(OUTFILE + ".nhaplo", "w")
# for idx,SELECTION in enumerate(select_arr):
for disp_idx,DISPERSAL_DISTANCE in enumerate(disp_arr):

	if(SELECTION == 0):
		file_prefix = DIR + SUBDIR + "/" + str(DISPERSAL_DISTANCE) + "_sampling" + "/*"+ str(SAMPLE_SIZE) + "_" + SUFFIX + ".trees"
	else:
		file_prefix = DIR + SUBDIR + "/" + str(DISPERSAL_DISTANCE) + "_" + str(SELECTION) + "_" + str(DOMINANCE) + "_sampling" + "/*"+ str(SAMPLE_SIZE) + "_" + SUFFIX + ".trees"
	files = glob.glob(file_prefix)  # use files generated from sampling_strategies_inds.py
	files = np.random.choice(files, N)
	print(len(files)); print(file_prefix)
	gs_files = []
	nh_files = []
	
	for file in tqdm(files):
		ts = tskit.load(file)

		gs = []
		nh = []

		for idx, window in enumerate(window_sizes):
	
			start = genome_length/2 - window/2
			stop = genome_length/2 + window/2

			subts = ts.keep_intervals([[start,stop]])
			g = subts.genotype_matrix()
			k = [hash(g[:, i].tobytes()) for i in range(g.shape[1])]
			
			hc = np.array(sorted(collections.Counter(k).values(), reverse=True))
			
			gsindex = 1 - np.sum((hc/SAMPLE_SIZE)**n)
			nhaps = len(hc)
			gs.append(gsindex)
			nh.append(nhaps)

		gs_files.append(gs)
		nh_files.append(nh)
	
	gs_avg = np.mean(gs_files, axis = 0)
	gs_std = np.std(gs_files, axis = 0)

	nh_avg = np.mean(nh_files, axis = 0)
	nh_std = np.std(nh_files, axis = 0)

	for idx, window in enumerate(window_sizes):
	
		output_gs.write(str(DISPERSAL_DISTANCE) + "\t" + str(SELECTION) + "\t" + str(window) + "\t" + "{:.5f}".format(gs_avg[idx]) + "\t" + "{:.5f}".format(gs_std[idx]) + "\n")
		output_nh.write(str(DISPERSAL_DISTANCE) + "\t" + str(SELECTION) + "\t" + str(window) + "\t" + "{:.5f}".format(nh_avg[idx]) + "\t" + "{:.5f}".format(nh_std[idx]) + "\n")

output_nh.close()
output_gs.close()
