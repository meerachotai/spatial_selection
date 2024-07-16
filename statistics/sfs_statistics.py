#!/usr/bin/env python3

# PREFIX can be taken from .trees filename, sampled .trees file will be named similarly
# /home2/mnc42/sampling_strategies_inds_dest.py ${PREFIX} ${MU} ${SAMPLE_SIZE} "sampling" ${SEED} --tree -mid
# /home2/mnc42/sfs.py /home2/mnc42/ spatial_WF_grid_sweep_from_neutral /home2/mnc42/spatial_WF_grid_sweep_from_neutral/sfs_window -rd

import argparse
import sys
import os
import numpy as np
import tskit
import msprime
import pyslim
import glob
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("selection", help="selection coef / 0 if neutral", type = float)
parser.add_argument("window", help="window size", type = float)
parser.add_argument("dir", help="main directory", type = str)
parser.add_argument("subdir", help="subdirectory", type = str)
parser.add_argument("out", help="outfile (including path)", type = str, default = "out")
parser.add_argument("nfiles", help="number of files", type = int)
parser.add_argument("-m","--mid", help="midpoint chromosome sampling", action = "store_true")
parser.add_argument("-rd","--random_diploid", help="random individual sampling", action = "store_true")

args = parser.parse_args()
SELECTION = args.selection
window_size = args.window
DIR=args.dir
SUBDIR=args.subdir
OUT=args.out
N=args.nfiles

# DIR = "/home2/mnc42/"
# SUBDIR = "spatial_WF_grid_sweep_from_neutral"
# SUFFIX = "_sampling_rh"

if(args.mid):
	SUFFIX = "sampling_mid"
if(args.random_diploid):
	SUFFIX = "sampling_rd"	

disp_arr = [0.015, 0.02, 0.04, 0.08, 0.1, 0.5]

# N = 1000 # number of replicates

# SFS and Tajima's D computation for selective sweeps

output = open(OUT + "_" + str(SELECTION) + "_" + str(window_size) + "_" + SUFFIX + ".sfs", "w")
output_d = open(OUT + "_" + str(SELECTION) + "_" + str(window_size) + "_" + SUFFIX + ".tajimad", "w")

START = (genome_length / 2) - (window_size/2)
STOP = (genome_length / 2) + (window_size/2)

if(window_size != genome_length):
	windows = [0, START, STOP, genome_length]
else:
	windows = [0, genome_length]
	
for DISPERSAL_DISTANCE in disp_arr:

	if(SELECTION == 0):
		file_prefix = DIR + SUBDIR + "/" + str(DISPERSAL_DISTANCE) + "_sampling" + "/*"+ str(SAMPLE_SIZE) + "_" + SUFFIX + ".trees"
	else:
		file_prefix = DIR + SUBDIR + "/" + str(DISPERSAL_DISTANCE) + "_" + str(SELECTION) + "_" + str(DOMINANCE) + "_sampling" + "/*"+ str(SAMPLE_SIZE) + "_" + SUFFIX + ".trees"

	print(file_prefix)
	files = glob.glob(file_prefix)
	print(len(files))
	files = np.random.choice(files, N)

	sfsfiles = []
	d = []

	for file in tqdm(files):
		ts = tskit.load(file)
		sfs = ts.allele_frequency_spectrum(windows = windows, polarised=True, span_normalise=False, mode = "site")
		seg = ts.segregating_sites(windows = windows, span_normalise = False)
		td = ts.Tajimas_D(windows = windows, mode = "site")

		if(window_size != 1e7):
			sfs = sfs[1]; seg = seg[1]; td = td[1]
		else:
			sfs = sfs[0]; seg = seg[0]; td = td[0]
		sfs = sfs / seg
		sfsfiles.append(sfs)
		d.append(td)

	sfs_mean = np.mean(sfsfiles,axis = 0)
	sfs_std = np.std(sfsfiles,axis = 0)

	d_avg = np.mean(d)
	d_std = np.std(d)

	output.write("AVG\t" + str(DISPERSAL_DISTANCE) + "\t" + str(SELECTION) + "\t" + '\t'.join("{:.10f}".format(item) for item in sfs_mean) + "\t" + str(len(sfsfiles)) + "\n")
	output.write("STD\t" + str(DISPERSAL_DISTANCE) + "\t" + str(SELECTION) + "\t" + '\t'.join("{:.10f}".format(item) for item in sfs_std) + "\t" + str(len(sfsfiles)) + "\n")

	output_d.write(str(DISPERSAL_DISTANCE) + "\t" + str(SELECTION) + "\t" + str(window_size) + "\t" + "{:.10f}".format(d_avg) + "\t" + "{:.10f}".format(d_std) + "\t" + str(len(files)) + "\n")

output.close()
output_d.close()
    