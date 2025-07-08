#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import random 

parser = argparse.ArgumentParser()

parser.add_argument("inprefix", help="inprefix (with path)", type = str)
parser.add_argument("replicate", help="replicate", type = int)
parser.add_argument("d", help="dispersal", type = float)

args = parser.parse_args()

inprefix = args.inprefix
DISPERSAL = args.d
REPLICATE = args.replicate

prefix = inprefix + "_" + str(REPLICATE)
file = prefix + "_2_" + str(DISPERSAL) + "_parentdist.txt"
print("reading", file)
d = pd.read_csv(file, sep = " ", header = None, index_col = None, names = ["offspring", "offspringx","offspringy","parent1","parent1x","parent1y","parent2","parent2x","parent2y"])

file = prefix + "_3_" + str(DISPERSAL) + "_fecundity.txt"
print("reading", file)
df = pd.read_csv(file, sep = " ", header = None, index_col = None, names = ["offspring", "offspringx","offspringy","fecundity"])

merged_d = d.merge(df, on=['offspringx', 'offspringy'])

# p1_d = (merged_d["offspringx"] - merged_d["parent1x"])**2
# p2_d = (merged_d["offspringx"] - merged_d["parent2x"])**2

# same as below:
p1_d_edges_x = abs((1 - (merged_d["offspringx"])) - merged_d["parent1x"])
p1_d_reg_x = abs(merged_d["offspringx"] - merged_d["parent1x"])
p1_d_x = np.minimum(p1_d_edges_x, p1_d_reg_x)

p1_d_edges_y = abs((1 - (merged_d["offspringy"])) - merged_d["parent1y"])
p1_d_reg_y = abs(merged_d["offspringy"] - merged_d["parent1y"])
p1_d_y = np.minimum(p1_d_edges_y, p1_d_reg_y)

euclidean_dist_p1 = np.sqrt((p1_d_x ** 2) + (p1_d_y ** 2))

p2_d_edges_x = abs((1 - (merged_d["offspringx"])) - merged_d["parent2x"])
p2_d_reg_x = abs(merged_d["offspringx"] - merged_d["parent2x"])
p2_d_x = np.minimum(p2_d_edges_x, p2_d_reg_x)

p2_d_edges_y = abs((1 - (merged_d["offspringy"])) - merged_d["parent2y"])
p2_d_reg_y = abs(merged_d["offspringy"] - merged_d["parent2y"])
p2_d_y = np.minimum(p2_d_edges_y, p2_d_reg_y)

euclidean_dist_p2 = np.sqrt((p2_d_x ** 2) + (p2_d_y ** 2))

average_parental_euclidean_d = (euclidean_dist_p1 + euclidean_dist_p2) / 2

mean_average_parental_euclidean_d = np.mean(average_parental_euclidean_d)

write_file = inprefix + "_" + str(DISPERSAL) + "_euclidean_dist.txt"

print("writing in file", write_file)
with open(write_file, "a") as effd_file:
    effd_file.write(str(REPLICATE) + "," + f"{mean_average_parental_euclidean_d:.20f}\n")

# write_file = inprefix + "_" + str(DISPERSAL) + "_averaged_v2.txt"
# print("writing in file", write_file)
# sampled_average_d = random.sample(list(average_d), 100)
# 
# with open(write_file, "a") as file:
#     file.write(str(REPLICATE) + ",")
#     for number in sampled_average_d:
#         file.write(f"{number:.20f},")
#     file.write("\n")