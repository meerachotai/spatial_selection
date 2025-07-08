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
matplotlib.rcParams['figure.dpi'] = 300

import pandas as pd

prefix = "data/supplement_data/distance/spatial_WF_grid_neutral"
eucd = []; eucd_std = []
num_inds = 10000

disp_array = [0.015, 0.02, 0.04, 0.1]
for idx,DISPERSAL in enumerate(disp_array):
    file = prefix + "_" + str(DISPERSAL) + "_euclidean_dist.txt"
    d = pd.read_csv(file, sep = ",", header = None, index_col = 0)
    eucd.append(d.iloc[:,0].mean())
    eucd_std.append(d.iloc[:,0].std())

fig, ax = plt.subplots(figsize = [2,2])

ax.plot(disp_array, eucd, marker = "o", color = "darkblue", linestyle = "dotted")
ax.errorbar(disp_array, eucd, yerr=eucd_std, capsize=1, color = "darkblue", linestyle = "dotted")

lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]

x = np.array([0] + disp_array + [0.12])
# now plot both limits against eachother
ax.plot(x, x/2, '--', color = "black", alpha=0.75, zorder=0, label = "y = x/2")

ax.set_xticks(disp_array)
ax.set_yticks(np.arange(0,0.09,0.02))
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.setp(ax.xaxis.get_majorticklabels(), rotation=55, ha = "right",rotation_mode="anchor", fontsize = 4.8)

ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.setp(ax.yaxis.get_majorticklabels(), rotation=0, ha = "right",rotation_mode="anchor", fontsize = 4.8)

ax.set_xlabel(r"dispersal ($d$)")
ax.set_ylabel(r"average euclidean distance")

ax.set_xlim([0,0.12])
ax.set_ylim([0,0.08])

handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor=(0.4, 0.85))

fig.savefig("figS1.pdf", bbox_inches = "tight", transparent = True)

