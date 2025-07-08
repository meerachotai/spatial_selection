#!/usr/bin/env python3

import tskit
import numpy as np
import glob
import pyslim
from tqdm import tqdm 
import pandas as pd
import argparse 
import random 
from sklearn.decomposition import PCA
import msprime

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
MU = 1e-8
REPLICATES = 1
min_val = 0; max_val = 0.1
sigma_arr=[0.015, 0.02, 0.04, 0.08, 0.1, 0.5]

for disp_idx,DISPERSAL_DISTANCE in enumerate(sigma_arr):

	for rep in range(1,REPLICATES + 1):
		file_prefix = DIR + "/" + SUBDIR + "/" + str(DISPERSAL_DISTANCE) + "/" + str(DISPERSAL_DISTANCE) + "_" + str(rep) + "*"+ str(MAXGEN) + ".trees"
		file = glob.glob(file_prefix)[0]
	
		print(file)
		nts = tskit.load(file)
		
# 		inds = np.random.choice(nts.num_individuals, int(SAMPLE_SIZE/2), replace = False)
# 		samps = []
# 		for i in inds:
# 			samps.extend(nts.individual(i).nodes) # using diploid individuals

		x = np.array([nts.individual(i).location[0] for i in range(nts.num_individuals)])
		y = np.array([nts.individual(i).location[1] for i in range(nts.num_individuals)])
		
		y_choose = list(np.nonzero((y > 0.45) & (y < 0.55))[0])
		
		x_choose_1 = list(np.nonzero((x > 0.2) & (x < 0.3))[0])
		inds_to_choose_1 = list(set(x_choose_1) & set(y_choose))
		inds_1 = np.random.choice(inds_to_choose_1, int(SAMPLE_SIZE/2), replace = False)
		
		x_choose_2 = list(np.nonzero((x > 0.7) & (x < 0.8))[0])
		
		inds_to_choose_2 = list(set(x_choose_2) & set(y_choose))
		inds_2 = np.random.choice(inds_to_choose_2, int(SAMPLE_SIZE/2), replace = False)
		
		samps = []
		loc_x = []; loc_y = []
		for i in inds_1:
			samps.extend(nts.individual(i).nodes) # using diploid individuals
			loc_x.append(nts.individual(i).location[0])
			loc_y.append(nts.individual(i).location[1])
		
		for i in inds_2:
			samps.extend(nts.individual(i).nodes) # using diploid individuals
			loc_x.append(nts.individual(i).location[0])
			loc_y.append(nts.individual(i).location[1])
		
		subsample_nodes=np.sort(np.array(samps))
		o=nts.simplify(subsample_nodes)
		
		next_id = pyslim.next_slim_mutation_id(nts)
		ts = msprime.sim_mutations(o,rate=MU,model=msprime.SLiMMutationModel(type=0, next_id=next_id),keep=True)
		
		geno = ts.genotype_matrix().T
		print(geno.shape)
		pca = PCA(n_components = 2)
		pca_fit = pca.fit_transform(geno)
		
		pca_df = pd.DataFrame(data = pca_fit, columns = ['PC1', 'PC2'])
		print(pca_fit.shape)
		pca_df.insert(2, "label",np.concatenate((np.zeros(SAMPLE_SIZE), np.ones(SAMPLE_SIZE))),allow_duplicates=False)
		
		pca_df.to_csv(OUT + "_" + str(DISPERSAL_DISTANCE) + "_" + str(rep) + ".pca", mode='w')
		
		loc_df = pd.DataFrame(data = np.array([loc_x, loc_y]).T, columns = ['x', 'y'])
		
		loc_df.insert(2, "label",np.concatenate((np.zeros(int(SAMPLE_SIZE/2)), np.ones(int(SAMPLE_SIZE/2)))),allow_duplicates=False)
		loc_df.to_csv(OUT + "_" + str(DISPERSAL_DISTANCE) + "_" + str(rep) + ".loc", mode='w')
		
# 		break
	