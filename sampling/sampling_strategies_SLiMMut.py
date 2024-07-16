#!/usr/bin/env python3

import sys
import os
import numpy as np
import tskit
import msprime
import pyslim
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("source_prefix", help="prefix of file (soure)", type = str)
parser.add_argument("dest_prefix", help="prefix of file (dest)", type = str)
parser.add_argument("mu", help="mutation rate", type = float)
parser.add_argument("n", help="sample size", type = int)
parser.add_argument("out", help="outsuffix to add onto prefix for output files", type = str)
parser.add_argument("seed", help="", type = int)
parser.add_argument("--vcf", help="output vcf", action = "store_true")
parser.add_argument("--tree", help="output tree", action = "store_true")
parser.add_argument("-m","--mid", help="midpoint chromosome sampling", action = "store_true")
parser.add_argument("-rd","--random_diploid", help="random individual sampling", action = "store_true")

args = parser.parse_args()
SPREFIX=args.source_prefix
DPREFIX=args.dest_prefix
MU=args.mu 
SAMPLE_SIZE=args.n
OUTSUFFIX=args.out
SEED=args.seed

def midpoint(nts,x,y, SAMPLE_SIZE, SEED):
	print("performing mid-point sampling")
	x_choose = list(np.nonzero((x > 0.45) & (x < 0.55))[0])
	y_choose = list(np.nonzero((y > 0.45) & (y < 0.55))[0])
	inds_to_choose = list(set(x_choose) & set(y_choose))

	np.random.seed(seed=SEED)
	inds = np.random.choice(inds_to_choose, int(SAMPLE_SIZE/2), replace = False) # outputs individual indices
	
	samps = []
	for i in inds:
		samps.extend(nts.individual(i).nodes) # using diploid individuals
	
	outfile = DPREFIX + "_" + str(SAMPLE_SIZE) + "_" + OUTSUFFIX + "_mid"
	
	return samps, outfile

def random_diploid(nts, SAMPLE_SIZE, SEED):
	print("performing random-diploid sampling")
	np.random.seed(seed=SEED)
	inds = np.random.choice(nts.num_individuals, int(SAMPLE_SIZE/2), replace = False) # outputs individual indices

	samps = []
	for i in inds:
		samps.extend(nts.individual(i).nodes) # using diploid individuals

	inds_out = [nts.node(i).individual for i in samps]
	assert set(inds) == set(inds_out)
	
	outfile = DPREFIX + "_" + str(SAMPLE_SIZE) + "_" + OUTSUFFIX + "_rd"
	
	return samps, outfile


def write_output(samps, outfile):
	subsample_nodes=np.sort(np.array(samps))
	o=nts.simplify(subsample_nodes)

	next_id = pyslim.next_slim_mutation_id(nts)
	ts = msprime.sim_mutations(o,rate=MU,model=msprime.SLiMMutationModel(type=0, next_id=next_id),keep=True)

# 	ts = msprime.sim_mutations(o, rate=MU, random_seed=SEED, keep = True)
	
	if(args.tree):
		print("writing trees file:", outfile)
		ts.dump(outfile + ".trees")

	if(args.vcf):
		print("writing vcf file:", outfile)
		vcfts = pyslim.generate_nucleotides(ts)
		vcfts = pyslim.convert_alleles(vcfts)

		inds = np.unique([nts.node(i).individual for i in samps])
		indv_names = [f"tsk_{i}indv" for i in range(len(inds))]
		with open(outfile + ".vcf", "w") as vcf_file:
			vcfts.write_vcf(vcf_file, individual_names=indv_names, isolated_as_missing = False)



file = SPREFIX + ".trees"
print(SPREFIX, MU, SAMPLE_SIZE)

nts = tskit.load(file)
x = np.array([nts.individual(i).location[0] for i in range(nts.num_individuals)])
y = np.array([nts.individual(i).location[1] for i in range(nts.num_individuals)])

if(args.mid):
	samps, outfile = midpoint(nts,x,y, SAMPLE_SIZE, SEED)
	print("number of samples =", len(samps))
	write_output(samps, outfile)
	
if(args.random_diploid):
	samps, outfile = random_diploid(nts,SAMPLE_SIZE, SEED)
	write_output(samps, outfile)


# next_id = pyslim.next_slim_mutation_id(nts)
# ts = msprime.sim_mutations(o,rate=MU,model=msprime.SLiMMutationModel(type=0, next_id=next_id),keep=True)
# nts = pyslim.generate_nucleotides(ts)
# nts = pyslim.convert_alleles(nts)
# indv_names = [f"tsk_{i}indv" for i in range(SAMPLE_SIZE)]
# with open(PREFIX + ".vcf", "w") as vcf_file:
#     nts.write_vcf(vcf_file, individual_names=indv_names)