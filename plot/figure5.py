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

fig, axs = plt.subplots(nrows=2, ncols = 2, figsize = [3.2,4])
gs = axs[1, 0].get_gridspec()
# remove the underlying axes
for ax in axs[1, :]:
    ax.remove()
axbig = fig.add_subplot(gs[1, :])

# -------------------------------Figure 5A-------------------------------------------------------------
pi_neutral = [0.000626,0.000407]

N = 10000
r = 1e-8
g = int(1e7)
center = int(g/2)

def exp_pi(SELECTION, DOMINANCE):
    calc_x = [abs(i - center) for i in range(g)] # distance from center
    x = [i for i in range(g)]
    sh = SELECTION * DOMINANCE # heterozygous s
    ex = -(2 * r * np.array(calc_x)) / sh
    y = 1 - ((4 * N * sh) ** ex)
    return x, y
    
# getting a subset
# output from sampled_pi_windows.py
import warnings

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

DIR = "data/"
SUBDIR = "sweep_size"

DISPERSAL_ARRAY = [0.015, 0.5]
SELECTION_ARRAY = [0.01, 0.1]
DOMINANCE = 0.5
SUFFIX = "sampling_rd"
average_interval = 10

colors = sns.diverging_palette(220, 20, s=150, l=45, n=2)
for axes_idx, SELECTION in enumerate(SELECTION_ARRAY):
    for idx_disp, DISPERSAL_DISTANCE in enumerate(DISPERSAL_ARRAY):
        idd = str(DISPERSAL_DISTANCE) + "_" + str(SELECTION) + "_" + str(DOMINANCE) 
        file_prefix =  DIR + SUBDIR + "/" + idd
        if SELECTION == 0.01:
            d = pd.read_csv(DIR + SUBDIR + "/" + "window_1e6_500_" + idd + "_" + SUFFIX + ".pi", delim_whitespace = True, header = 0)
        else:
            d = pd.read_csv(DIR + SUBDIR + "/" + "window_" + idd + "_" + SUFFIX + ".pi", delim_whitespace = True, header = 0)
        averages_d = pd.DataFrame()
        for i in range(0, d.shape[1], average_interval):
            subset = d.iloc[:, i:i+average_interval]
            averages = subset.mean(axis=1)
            averages_d[d.columns[i]] = averages
        mean = averages_d.mean(axis = 0)#/window_length # per bp for 10kbp window
        windows = averages_d.columns.astype(float)
        mean_adj = mean/pi_neutral[idx_disp]
        
        if(DISPERSAL_DISTANCE == 0.5):
            DISPERSAL_DISTANCE = 1
        LABEL = r"$d$ = " + "{:.3f}".format(DISPERSAL_DISTANCE)
        axs[0,axes_idx].scatter(windows/1000, mean_adj,marker='o', color = colors[idx_disp], label = LABEL, alpha = 0.5)
        axs[0,axes_idx].set_title(r"$s$ = " + "{:.2f}".format(SELECTION))
    axs[0,axes_idx].set_xlabel("position (in kbp)")
    x, y = exp_pi(SELECTION, DOMINANCE)
    axs[0,axes_idx].plot(np.array(x[(int(windows[0])):(int(windows[-1]))])/1000, np.array(y[(int(windows[0])):(int(windows[-1]))]), color = "grey", linestyle = "--")

axs[0,0].set_ylim([0,1.21])
axs[0,1].set_ylim([0,1.21])

axs[0,0].set_ylabel(r"$\frac{\pi(x)}{\pi_0}$")
axs[0,1].get_yaxis().set_visible(False)
# https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.legend.html
axs[0,0].legend(bbox_to_anchor=(0.48, 0.23), loc='upper right',markerscale=0.7, handletextpad = 0.03, borderpad = 0.2)

# -------------------------------Figure 5B-------------------------------------------------------------
colors = sns.diverging_palette(220, 20, s=150, l=45, n=6)

DIR = "data/"
SUBDIR = "roh"
GENERATION = 500000
LENGTH = "1e7"
axes_idx = 0

sigma_arr=[0.015, 0.02, 0.04, 0.08, 0.1, 0.5]
SAMPLE_SIZE = 100
REPLICATES = 3
SNP_KB = "_20_200" # CHANGE THIS AS NEEDED
SUFFIX = str(GENERATION) + "_100_sampling_rd" + SNP_KB

bins = np.arange(0, 11000, 200)

sigma_arr = [0.015, 0.02, 0.04, 0.1,0.5]
for idx,DISPERSAL_DISTANCE in enumerate(sigma_arr):
    prefix =  DIR + SUBDIR + "/" + str(DISPERSAL_DISTANCE) + "/" + str(DISPERSAL_DISTANCE)
    print(prefix)
    files = glob.glob(prefix + "*" + "_" + SUFFIX + ".hom")
    
    yall = []
    for file in files:
        print(file)
        d = pd.read_csv(file, delim_whitespace = True, header = 0)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        d['window'] = pd.cut(d['KB'], bins = bins, labels = bin_centers)
        print(np.sum(d['KB']))
        
        y = []
        for x in bin_centers:
            subd = d[d['window'] == x]
            y.append((np.sum(subd['KB'])*1e3 * 100)/(1*1e7*50)) # divided by length of genomes 50 individuals x 10Mb
        yall.append(y)
    ymean = np.mean(yall, axis = 0) # average across replicates
    print(sum(ymean))
    ycs = np.cumsum(ymean)
    print("max:",max(ycs))
    
    if(DISPERSAL_DISTANCE == 0.5):
        DISPERSAL_DISTANCE = 1
    axbig.plot(bin_centers/1000, ycs, linestyle='-', color=colors[idx], label = r"$d$ = " + "{:.3f}".format(DISPERSAL_DISTANCE))

axbig.set_xlabel("length of ROH (Mb)")
axbig.set_ylabel("% of the genome")
axbig.set_xlim([0,10]) # 10
axbig.set_ylim([0,25]) # 10
# axbig.ticklabel_format(axis='x', style='sci', scilimits=(3,3), useMathText = True)


axbig.legend(bbox_to_anchor=(0.34, 1), loc='upper right', borderpad = 0.25)

ax = axs.flatten()
for i, label in enumerate(('A', '')):
    ax[i].text(-0.35, 1.12, label, transform=ax[i].transAxes, fontsize = 8,va='top', ha='right', weight = "bold")

axbig.text(-0.17, 1.12, "B", transform=axbig.transAxes, fontsize = 8,va='top', ha='right', weight = "bold")


# -----------------------------------------------------------------------
fig.tight_layout()
fig.savefig("fig5" + ".pdf",transparent=True, bbox_inches = "tight")#, pad_inches = 0)  

