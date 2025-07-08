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
matplotlib.rcParams['lines.linewidth'] = 1
matplotlib.rcParams['lines.markersize'] = 2

from matplotlib.lines import Line2D

fig, axs = plt.subplots(nrows=1, ncols = 1, figsize = [2.5,2])
# gs = axs[1, 0].get_gridspec()
# remove the underlying axes
# for ax in axs[1, :]:
    # ax.remove()
# axbig = fig.add_subplot(gs[1, :])

# -----------------------------------------------------------------------
pi_neutral = [61385.148280,39910.496561,39910.496561]

N = 10000
r = 1e-8
g = int(1e7)
center = int(g/2)

def exp_pi(SELECTION, DOMINANCE,x, center):
    calc_x = [abs(i - center) for i in x] # distance from center
    sh = SELECTION * DOMINANCE # heterozygous s
    ex = -(2 * r * np.array(calc_x)) / sh
    y = 1 - ((4 * N * sh) ** ex)
    return y
    
# getting a subset
# output from sampled_pi_windows.py
import warnings

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

DIR = "data/supplement_data"
SUBDIR = "sweep_center"

DISPERSAL_ARRAY = [0.015, 0.5]
SELECTION = 0.1
DOMINANCE = 0.5
SUFFIX = "sampling_rd"
average_interval = 20
slide_interval = 10

colors = sns.diverging_palette(220, 20, s=150, l=45, n=2)
colors = colors + ["coral"]

for idx_disp, DISPERSAL_DISTANCE in enumerate(DISPERSAL_ARRAY):
    idd = str(DISPERSAL_DISTANCE) + "_" + str(SELECTION) + "_" + str(DOMINANCE) 
    file_prefix =  DIR + "/" + SUBDIR + "/" + idd

    d = pd.read_csv(DIR + "/" + SUBDIR + "/" + "window_center_" + idd + "_" + SUFFIX + ".branch_pi", sep = "\t", header = 0)
    
    pd_cols = []
    cols = np.arange(0, d.shape[1]-(average_interval-slide_interval)-1, average_interval-slide_interval)
    for i in cols:
        subset = d.iloc[:, i:i+average_interval]
        # print(subset.shape)
        averages = subset.mean(axis=1)
        pd_cols.append(averages)
        if(i < cols[3]):
            print(d.columns[i], d.columns[i+average_interval])
    averages_d = pd.concat(pd_cols, axis = 1)
    averages_d.columns = d.columns[cols]
    # print(d.columns[cols][:4])
    
    mean = averages_d.mean(axis = 0)#/window_length # per bp for 10kbp window
    windows = averages_d.columns.astype(float)

    mean_adj = mean/pi_neutral[idx_disp]
    
    if(DISPERSAL_DISTANCE == 0.5):
        DISPERSAL_DISTANCE = 1
    LABEL = r"$d$ = " + "{:.3f}".format(DISPERSAL_DISTANCE)
    axs.scatter((windows - 5e6)/1000, mean_adj,marker='o', color = colors[idx_disp], label = LABEL, alpha = 0.5)

# axs.set_title(r"$s$ = " + "{:.2f}".format(SELECTION))
axs.set_xlabel("position relative to adaptive allele (in kbp)")

# ---- add SGV ------------
FREQUENCY = 0.05
idx_disp = 2
idd = str(FREQUENCY) + "_" + str(SELECTION) + "_" + str(DOMINANCE) 
file_prefix =  DIR +"/" + SUBDIR + "/" + idd

d = pd.read_csv(DIR +"/" + SUBDIR + "/" + "window_center_" + idd + "_" + SUFFIX + ".branch_pi", sep = "\t", header = 0)

pd_cols = []
cols = np.arange(0, d.shape[1]-(average_interval-slide_interval)-1, average_interval-slide_interval)
for i in cols:
    subset = d.iloc[:, i:i+average_interval]
    # print(subset.shape)
    averages = subset.mean(axis=1)
    pd_cols.append(averages)
    # if(i < cols[3]):
        # print(d.columns[i], d.columns[i+average_interval])
averages_d = pd.concat(pd_cols, axis = 1)
averages_d.columns = d.columns[cols]
# print(d.columns[cols][:4])

mean = averages_d.mean(axis = 0)#/window_length # per bp for 10kbp window
windows = averages_d.columns.astype(float)
mean_adj = mean/pi_neutral[idx_disp]

LABEL = r"$x_0$ = " + "{:.3f}".format(FREQUENCY)
axs.scatter((windows - 5e6)/1000, mean_adj,marker='o', color = colors[idx_disp], label = LABEL, alpha = 0.5)


# ------ plot expected --------
x = [float(i) for i in d.columns]
print(x[:3])
y = exp_pi(SELECTION, DOMINANCE, x,center)
print(y[:3])
cols = np.arange(0, d.shape[1]-(average_interval-slide_interval)-1, average_interval-slide_interval)
sliding_y = []
for i in cols:
    sliding_y.append(y[i:(i+average_interval)].mean())
axs.plot((windows - 5e6)/1000, sliding_y, color = "grey", linestyle = "--")

axs.set_ylabel(r"mean tree height (relative to neutral)")
# https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.legend.html
axs.legend(bbox_to_anchor=(0.33, 0.9), loc='upper right',markerscale=0.7, handletextpad = 0.03, borderpad = 0.2)

axs.set_ylim(0,1)
axs.set_xticks(np.arange(-50,50+5,10))

fig.tight_layout()
fig.savefig("figS7.pdf",transparent=True, bbox_inches = "tight")#, pad_inches = 0)  
