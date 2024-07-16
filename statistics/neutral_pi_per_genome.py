#!/usr/bin/env python3

import tskit
import numpy as np
import glob
import pyslim
import seaborn as sns
from tqdm import tqdm 
import pandas as pd
import argparse 
import random 

parser = argparse.ArgumentParser()
parser.add_argument("dir", help="main directory", type = str)
parser.add_argument("subdir", help="subdirectory", type = str)
parser.add_argument("out", help="outfile (including path)", type = str, default = "out")
parser.add_argument("-m","--mid", help="midpoint chromosome sampling", action = "store_true")
parser.add_argument("-rd","--random_diploid", help="random individual sampling", action = "store_true")

args = parser.parse_args()
DIR=args.dir
SUBDIR=args.subdir
OUT=args.out

if(args.mid):
	SUFFIX = "sampling_mid"
if(args.random_diploid):
	SUFFIX = "sampling_rd"	

SAMPLE_SIZE = 100
N = 100

sigma_arr=[0.015, 0.02, 0.04, 0.08, 0.1, 0.5]

output = open(OUT + "_" + SUFFIX + ".genome_pi", "w")

INTERVAL = 1e5
windows = np.arange(0,1e7+INTERVAL, INTERVAL)
window_idx = np.random.randint(0, len(windows)-1)

for disp_idx,DISPERSAL_DISTANCE in enumerate(sigma_arr):
	file_prefix = DIR + "/" + SUBDIR + "/" + str(DISPERSAL_DISTANCE) + "_sampling" + "/*"+ str(SAMPLE_SIZE) + "_" + SUFFIX + ".trees"
	print(file_prefix)
	files = glob.glob(file_prefix)
	random.shuffle(files)
	p_list = []
	for file in tqdm(files[:N]):
		ts = tskit.load(file)
		inds = []
		for ind in ts.individuals():
			nodes = ind.nodes
			inds.append(nodes.tolist())
		diversity = ts.diversity(sample_sets = inds,windows = windows,span_normalise = False)
		diversity = diversity[window_idx]/INTERVAL
		p_list.extend(diversity)
	output.write(str(DISPERSAL_DISTANCE) + "\t" + '\t'.join("{:.10f}".format(item) for item in p_list) + "\n")

output.close()
