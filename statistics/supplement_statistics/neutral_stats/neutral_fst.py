#!/usr/bin/env python3

import tskit
import numpy as np
import glob
import pyslim
from tqdm import tqdm 
import pandas as pd
import argparse 
import random 

parser = argparse.ArgumentParser()
parser.add_argument("dir", help="main directory", type = str)
parser.add_argument("subdir", help="subdirectory", type = str)
parser.add_argument("out", help="outfile (including path)", type = str, default = "out")

args = parser.parse_args()
DIR=args.dir
SUBDIR=args.subdir
OUT=args.out

MAXGEN = 500000

SAMPLE_SIZE = 100
min_val = 0; max_val = 0.1
sigma_arr=[0.015, 0.02, 0.04, 0.08, 0.1, 0.5]

output = open(OUT + ".fst", "w")


for disp_idx,DISPERSAL_DISTANCE in enumerate(sigma_arr):
	file_prefix = DIR + "/" + SUBDIR + "/" + str(DISPERSAL_DISTANCE) + "/" + str(DISPERSAL_DISTANCE) + "*"+ str(MAXGEN) + ".trees"
	print(file_prefix)
	files = glob.glob(file_prefix)[:3]
	fst = []
	for file in tqdm(files):
		nts = tskit.load(file)
		x = np.array([nts.individual(i).location[0] for i in range(nts.num_individuals)])
		y = np.array([nts.individual(i).location[1] for i in range(nts.num_individuals)])
		
		y_choose = list(np.nonzero((y > 0.45) & (y < 0.55))[0])
		
		x_choose_1 = list(np.nonzero((x > 0.2) & (x < 0.3))[0])
		inds_to_choose_1 = list(set(x_choose_1) & set(y_choose))
		inds_1 = np.random.choice(inds_to_choose_1, int(SAMPLE_SIZE/2), replace = False)
		samps_1 = []
		for i in inds_1:
			samps_1.extend(nts.individual(i).nodes) # using diploid individuals
		 
		x_choose_2 = list(np.nonzero((x > 0.7) & (x < 0.8))[0])
		
		inds_to_choose_2 = list(set(x_choose_2) & set(y_choose))
		inds_2 = np.random.choice(inds_to_choose_2, int(SAMPLE_SIZE/2), replace = False)
		
		samps_2 = []
		for i in inds_2:
			samps_2.extend(nts.individual(i).nodes) # using diploid individuals
		
		fst.append(nts.Fst([samps_1, samps_2], mode = "branch"))

	output.write(str(DISPERSAL_DISTANCE) + "\t" + '\t'.join("{:.10f}".format(item) for item in fst) + "\n")
		
output.close()
