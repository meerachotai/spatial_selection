#!/usr/bin/env python3

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
parser.add_argument("dispersal", help="dispersal distance", type = float)
parser.add_argument("n", help="sample size", type = int)
parser.add_argument("dir", help="main directory", type = str)
parser.add_argument("subdir", help="subdirectory", type = str)
parser.add_argument("seed", help="seed", type = int)
parser.add_argument("out", help="outfile (including path)", type = str, default = "out")

args = parser.parse_args()
DISPERSAL_DISTANCE= args.dispersal
SAMPLE_SIZE = args.n
SEED = args.seed
DIR=args.dir
SUBDIR=args.subdir
OUT=args.out

genome_length = 1e7
generation = 500000
# if(DISPERSAL_DISTANCE == 0):
# 	DISPERSAL_DISTANCE = int(DISPERSAL_DISTANCE)

file_prefix = DIR + "/" + SUBDIR + "/" + str(DISPERSAL_DISTANCE)
print(file_prefix)
files = glob.glob(file_prefix + '/*' + str(generation) + '.trees')

print(files)
np.random.seed(SEED)
file = np.random.choice(files) # one replicate for now
nts = tskit.load(file)

output = open(OUT + "_" + str(DISPERSAL_DISTANCE) + ".tmrca", "w")
lengths = np.diff(nts.breakpoints(as_array = True)) / nts.sequence_length
print("number of trees =", nts.num_trees)
for idx,tree in tqdm(enumerate(nts.trees())):
    samps = [i for i in tree.samples() if nts.node(i).time == 0] # only present day samples
    np.random.seed(SEED)
    np.random.shuffle(samps)
    i = 0
    tmrca = []
    while(len(tmrca) < (SAMPLE_SIZE / 2)):
        mrca = tree.mrca(samps[i*2],samps[(i*2) + 1])
        if(mrca > 0): # if coalesced
            tmrca.append(tree.time(mrca))
        i = i + 1
    output.write(str(idx) + "\t" + str(int(lengths[idx] * genome_length)) + "\t" + '\t'.join("{:.10f}".format(item) for item in tmrca) + "\n")    
    