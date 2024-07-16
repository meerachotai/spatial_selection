#!/usr/bin/env python3

# PREFIX can be taken from .trees filename, sampled .trees file will be named similarly
# /home2/mnc42/pi_windows.py 0.05 /home2/mnc42/ spatial_WF_grid_sweep_from_neutral /home2/mnc42/spatial_WF_grid_sweep_from_neutral/window 2e6 2000 -rd
# /home2/mnc42/pi_windows.py 0.1 /home2/mnc42/ spatial_WF_grid_sweep_from_neutral /home2/mnc42/spatial_WF_grid_sweep_from_neutral/window 4e6 4000 -rd

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
parser.add_argument("selection", help="selection coeff", type = str)
parser.add_argument("dir", help="main directory", type = str)
parser.add_argument("subdir", help="subdirectory", type = str)
parser.add_argument("out", help="outfile (including path)", type = str, default = "out")
parser.add_argument("window", help="window size", type = float)
parser.add_argument("subwindow", help="subwindow size", type = float, default = 1e4)
parser.add_argument("nfiles", help="number of files", type = int)
parser.add_argument("-m","--mid", help="midpoint chromosome sampling", action = "store_true")
parser.add_argument("-rd","--random_diploid", help="random individual sampling", action = "store_true")


args = parser.parse_args()
SELECTION=args.selection
DIR=args.dir
SUBDIR=args.subdir
OUT=args.out
WINDOW_SIZE=args.window
SUB_WINDOW_SIZE=args.subwindow
N = args.nfiles # number of replicates

# DIR = "/home2/mnc42/"
# SUBDIR = "spatial_WF_grid_sweep_from_neutral"
# SUFFIX = "_sampling_rh"

if(args.mid):
	SUFFIX = "sampling_mid"
if(args.random_diploid):
	SUFFIX = "sampling_rd"	

DISPERSAL_ARRAY = [0.015, 0.5] # 0 is panmictic
DOMINANCE = 0.5

genome_length = 1e7

START = genome_length/2 - WINDOW_SIZE/2 
STOP =  genome_length/2 + WINDOW_SIZE/2

num_windows = int(WINDOW_SIZE / SUB_WINDOW_SIZE) # 10kb windows
windows = [0]
windows.extend(np.linspace(START, STOP, num_windows + 1))
windows.append(genome_length)

SAMPLE_SIZE = 100

piWindows = []

for DISPERSAL_DISTANCE in DISPERSAL_ARRAY:
	print(DISPERSAL_DISTANCE, SELECTION)
	file_prefix = DIR + SUBDIR + "/" + str(DISPERSAL_DISTANCE) + "_" + str(SELECTION) + "_" + str(DOMINANCE) + "_sampling"  + "/*"+ str(SAMPLE_SIZE) + "_" + SUFFIX + ".trees"
	files = glob.glob(file_prefix)
	files = np.random.choice(files, N)
	print(len(files))
	
	outfile = OUT + "_" + str(DISPERSAL_DISTANCE) + "_" + str(SELECTION) + "_" + str(DOMINANCE) + "_" + SUFFIX + ".pi"
	print(outfile)
	pi_output = open(outfile, "w")
	pi_output.write('\t'.join("{:.5f}".format(item) for item in windows[1:-1]) + "\n")
	pi = []
	for file in tqdm(files):
		nts = tskit.load(file)
		p = np.squeeze(nts.diversity(windows = windows, span_normalise = False))[1:-1] # ignore first, last windows
		p = p/SUB_WINDOW_SIZE
		pi_output.write('\t'.join("{:.5f}".format(item) for item in p) + "\n")
	pi_output.close()

