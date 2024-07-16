#!/usr/bin/env python3

import msprime
import pyslim
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("prefix", help="prefix of file", type = str)
parser.add_argument("mu", help="mutation rate", type = float)
# parser.add_argument("n", help="sample size", type = int)
parser.add_argument("seed", help="seed", type = int)
parser.add_argument("freq", help="starting frequency", type = float)
parser.add_argument("r", help="recombination", type = float)
parser.add_argument("l", help="genome length", type = float)
parser.add_argument("s", help="selection", type = float)
parser.add_argument("popsize", help="population size", type = int)

args = parser.parse_args()
OUTPREFIX=args.prefix
MU=args.mu
# SAMPLE_SIZE=args.n
SEED=args.seed
RECOMBINATION=args.r
SELECTION=args.s
POPSIZE=args.popsize
GENOME_LENGTH=args.l
STARTING_FREQUENCY=args.freq

WINDOW_SIZE = 500000

ts = msprime.sim_ancestry(
    samples=[msprime.SampleSet(POPSIZE)],
    population_size=POPSIZE,
    recombination_rate=RECOMBINATION,
    sequence_length=GENOME_LENGTH,
    random_seed=SEED,
)

# annotated tree sequence
ots = pyslim.annotate(ts, model_type="WF", tick=1, stage="late")

# this means the mutations will be of type “m2” in SLiM 
# (and, so you must initialize that mutation type in the recipe that loads this tree sequence in).

mut_model = msprime.SLiMMutationModel(type=2)
# annotated + mutated tree sequence
mts = msprime.sim_mutations(
            ots,rate=MU,
            model=mut_model,
            keep=True, 
            random_seed=SEED)

print(f"The tree sequence now has {mts.num_mutations} mutations, at "
      f"{mts.num_sites} distinct sites.")


geno = mts.genotype_matrix()
freq = geno.sum(axis = 1)

idx = np.where(freq == (STARTING_FREQUENCY * geno.shape[1]))[0]
pos_mut = []
for i in idx:
    if (len(mts.site(i).mutations) == 1):
        pos = mts.site(i).position
        start = pos - WINDOW_SIZE/2; stop = pos + WINDOW_SIZE/2
        if(stop < 1e7 or start > 0):
            pos_mut.append(mts.site(i).mutations[0].id)

# np.random.seed(seed=SEED)
# inds = np.random.choice(mts.num_individuals, int(SAMPLE_SIZE/2), replace = False)
# samps = []
# for i in inds:
#     samps.extend(mts.individual(i).nodes) # using diploid individuals
#     
# pos_mut = []
# for variant in mts.variants():
#     if(len(variant.site.mutations) == 1): # ignore stacked mutations for now
#         if((sum(variant.genotypes[samps])) == (STARTING_FREQUENCY * SAMPLE_SIZE)):
#             pos_mut.append(variant.site.mutations[0].id)
            
np.random.seed(SEED)
mut_id = np.random.choice(pos_mut)

tables = mts.dump_tables()
tables.mutations.clear()

m = mts.mutation(mut_id)
md_list = m.metadata["mutation_list"]
slim_ids = m.derived_state.split(",")
assert len(slim_ids) == len(md_list)
for sid, md in zip(slim_ids, md_list):
    md["selection_coeff"] = SELECTION
_ = tables.mutations.append(m.replace(metadata={"mutation_list": md_list}))

# tree sequence with only one selected mutation
mod_ts = tables.tree_sequence()

mod_ts.dump(OUTPREFIX + "_init.trees")

