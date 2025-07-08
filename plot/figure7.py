#!/usr/bin/env python3

import sys
import os
import numpy as np
import glob
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import seaborn as sns

font = {'size': 5.5}
matplotlib.rc('font', **font)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['lines.linewidth'] = 0.5
matplotlib.rcParams['lines.markersize'] = 2

from matplotlib.lines import Line2D

fig, axs1 = plt.subplots(nrows= 1, ncols = 1, figsize = [2,2])

# -------------------------------------------------------------------------------------------------------

DOMINANCE = 0.5
TOTAL_FILES = 10
SAMPLE_SIZE = 100
MU=1e-8
genome_length = 1e7
n = 2

window = 100000 # other window sizes plotted in supplement
SELECTION = 0.1

N = 1000; SAMPLE_SIZE = 100


colors = sns.diverging_palette(220, 20, s=150, l=45, n=2)

color_idx = 0

DIR = "data/"
SUBDIR = "haplotype_statistics"
PREFIX = "out"
SUFFIX = "sampling_rd"
file = DIR + SUBDIR + "/" + PREFIX + "_" + str(SELECTION) + "_" + str(N) + "_" + str(SAMPLE_SIZE) + "_" + SUFFIX + ".het"
print(file)
d = pd.read_csv(file, delim_whitespace = True, header = None,
               names = ["DISPERSAL", "SELECTION", "WINDOW", "HET_AVG", "HET_STD"])
d['DISPERSAL'] = d['DISPERSAL'].apply(lambda x: 1 if x == 0.5 else x)

LABEL = f"{window / 10**5:.1f}" + r"$\times 10^5$"
subset = d[(d["SELECTION"] == SELECTION) & (d["WINDOW"] == window)]
subset = subset.sort_values(by='DISPERSAL')
axs1.plot(subset["DISPERSAL"],subset["HET_AVG"], marker = "o", linestyle = "--", color = np.array(colors)[color_idx],label = LABEL)
axs1.errorbar(subset["DISPERSAL"], subset["HET_AVG"],yerr=subset["HET_STD"], capsize=1, color = np.array(colors)[color_idx], linestyle = "--")

# axs1.set_xlabel("hard sweep in spatial populations (starting frequency = 0.0001)\ndispersal",color=np.array(colors)[color_idx])
axs1.set_xlabel(r"dispersal ($d$)",color=np.array(colors)[color_idx])

axs1.set_xscale("log")
axs1.set_xticks(subset["DISPERSAL"])
axs1.tick_params(axis='x', rotation=55,labelcolor=np.array(colors)[color_idx])
axs1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

# -------------------------------------------------------------------------------------------------------

axs2 = axs1.twiny()
axs1.set_ylim([0.4,1])
axs1.set_ylabel("haplotype heterozygosity")
# -------------------------------------------------------------------------------------------------------

DIR = "data/"
SUBDIR = "haplotype_SGV_statistics"

color_idx = 1
x = np.arange(0, SAMPLE_SIZE + 1, 1)

SUFFIX = "sampling_rd"
file = DIR + SUBDIR + "/out_" + str(N) + "_" + str(SAMPLE_SIZE) + "_" + SUFFIX + "_" + "SGV_pop.het"
print(file)
d = pd.read_csv(file, delim_whitespace = True, header = None,
               names = ["FREQUENCY", "SELECTION", "WINDOW", "HET_AVG", "HET_STD"])
d = d[d['FREQUENCY'] != 0.04]

LABEL = f"{window / 10**5:.1f}" + r"$\times 10^5$"
subset = d[(d["SELECTION"] == SELECTION) & (d["WINDOW"] == window)]
subset = subset.sort_values(by='FREQUENCY')
axs2.plot(subset["FREQUENCY"],subset["HET_AVG"], marker = "o", color = np.array(colors)[color_idx],label = LABEL, linestyle = "--")
axs2.errorbar(subset["FREQUENCY"],subset["HET_AVG"], yerr=subset["HET_STD"], capsize=1, color = np.array(colors)[color_idx], linestyle = "--")

# axs2.set_xlabel("SGV soft sweeps in panmictic populations (dispersal = 1)\nstarting frequency",color=np.array(colors)[color_idx])
axs2.set_xlabel("starting frequency",color=np.array(colors)[color_idx])

axs2.set_xscale("log")
axs2.set_xticks(subset["FREQUENCY"])
axs2.tick_params(axis='x', rotation=75,labelcolor=np.array(colors)[color_idx])
axs2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

axs2.invert_xaxis()

fig.tight_layout()

fig.savefig("fig7.pdf", bbox_inches = "tight", transparent = True)