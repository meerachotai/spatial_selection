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
import seaborn as sns

gridspec = dict(hspace=0.1, width_ratios=[1, 1, 0.6, 1, 1], wspace = 0.05)

fig, axs = plt.subplots(nrows=2, ncols = 5, figsize = [5.6,3], gridspec_kw=gridspec)
for i in range(2):
    axs[i,2].set_visible(False)

#------------------------------------------------------------Figure S9A------------------------------------------------------------
DIR = "data/supplement_data/extend_sweeps"
PREFIX = "sfs_time"
ORIG_DIR = "data/sfs"

SAMPLE_SIZE = 100

SUFFIXES = ["sampling_rd","sampling_mid"]

colors = sns.diverging_palette(220, 20, s=150, l=45, n=4)
x = np.arange(0, SAMPLE_SIZE + 1, 1)
SELECTION = 0.1
dispersal_arr = [0.015,0.5]
time_arr = [0, 9000, 9900, 10000]
window_size = float(1000000)

show_time = 10000 - np.array(time_arr)

for sampling_idx, SUFFIX in enumerate(SUFFIXES):
    for disp_idx,DISPERSAL_DISTANCE in enumerate(dispersal_arr):
        file = DIR + "/" + PREFIX + "_" + str(SELECTION) + "_" + str(DISPERSAL_DISTANCE) + "_" + str(window_size) + "_" + SUFFIX + ".sfs"
        print(SUFFIX, DISPERSAL_DISTANCE, file)
        d = pd.read_csv(file, sep = "\t", header = None,
                   names = ["TYPE","TIME", "SELECTION"] + [i for i in range(SAMPLE_SIZE + 1)] + ["FILE_COUNT"])
        
        for color_idx,TIME in enumerate(time_arr):
            if(TIME == 10000):
                d0 = pd.read_csv(ORIG_DIR + "/sfs_window_" + str(SELECTION) + "_" + str(window_size) + "_" + SUFFIX + ".sfs", sep = "\t", header = None,
                       names = ["TYPE","DISPERSAL", "SELECTION"] + [i for i in range(SAMPLE_SIZE + 1)] + ["FILE_COUNT"])
                sfs = d0[(d0["DISPERSAL"] == DISPERSAL_DISTANCE) & (d0["TYPE"] == "AVG") & (d0["SELECTION"] == SELECTION)].iloc[:, 3:]
            else:
                sfs = d[(d["TIME"] == TIME) & (d["TYPE"] == "AVG") & (d["SELECTION"] == SELECTION)].iloc[:, 3:]
            label_time = show_time[color_idx]
            # LABEL = "$t$ = " + f"{label_time / 10**4:.2f}" + r"$\times 10^4$"
            LABEL = "$t$ = " + f"{label_time}"
            
            if(sampling_idx == 0):
                axs[sampling_idx,disp_idx].scatter(sfs.columns[1:-1]/100,sfs.iloc[0,1:-1].tolist(), color = np.array(colors)[color_idx], label = LABEL, s = 0.8,alpha = 0.8)
            else:
                axs[sampling_idx,disp_idx].scatter(sfs.columns[1:-1]/100,sfs.iloc[0,1:-1].tolist(), color = np.array(colors)[color_idx], s = 0.8, alpha = 0.8)

        # adding neutral
        dn = pd.read_csv(ORIG_DIR + "/sfs_window_0.0_" + str(float(1e7)) + "_" + SUFFIX + ".sfs", sep = "\t", header = None,
               names = ["TYPE","DISPERSAL", "SELECTION"] + [i for i in range(SAMPLE_SIZE + 1)] + ["FILE_COUNT"])
        sfsn = dn[(dn["DISPERSAL"] == DISPERSAL_DISTANCE) & (dn["TYPE"] == "AVG") & (dn["SELECTION"] == 0)].iloc[:, 3:]
        axs[sampling_idx,disp_idx].scatter(sfsn.columns[1:-1]/100,sfsn.iloc[0,1:-1].tolist(), color = "grey", s = 0.8,alpha = 0.5)#, label="neutral")
        # axs[sampling_idx,disp_idx].plot(sfsn.columns[1:-1]/100,sfsn.iloc[0,1:-1].tolist(), color = "grey", alpha = 0.5, label="neutral")

# -------------------------------------------------------------Figure S9B------------------------------------------------------------

DOMINANCE = 0.5
TOTAL_FILES = 10
SAMPLE_SIZE = 100
N = 1000

MU=1e-8
SEED = 5
genome_length = 1e7
n = 2

DIR = "data/supplement_data/extend_sweeps"
PREFIX = "out_time"
ORIG_DIR = "data/haplotype_statistics"

window_sizes = [10000, 30000,100000]

colors_all = sns.diverging_palette(145, 300, s=60, l = 50, n = len(SUFFIXES) * len(window_sizes))# + 2)
colors_local = colors_all[:3][::-1]
colors_global = colors_all[-3:]
colors = [colors_global] + [colors_local]

x = np.arange(0, SAMPLE_SIZE + 1, 1)
dispersal_arr = [0.015,0.5]
SELECTION = 0.1

for sampling_idx, SUFFIX in enumerate(SUFFIXES):
    for idx,DISPERSAL in enumerate(dispersal_arr):
        file = DIR + "/" + PREFIX + "_" + str(DISPERSAL) + "_" + str(SELECTION) + "_" + str(SAMPLE_SIZE) + "_" + SUFFIX + ".het"
        print(SUFFIX, DISPERSAL, file)#, os.path.getmtime(file))
        d = pd.read_csv(file, sep = "\t", header = None,
                   names = ["TIME", "SELECTION", "WINDOW", "HET_AVG", "HET_STD"])

        file0 = ORIG_DIR + "/out" + "_" + str(SELECTION) + "_" + str(N) + "_" + str(SAMPLE_SIZE) + "_" + SUFFIX + ".het"
        d0 = pd.read_csv(file0, sep = "\t", header = None,
                   names = ["DISPERSAL", "SELECTION", "WINDOW", "HET_AVG", "HET_STD"])
        
        # adding neutral
        filen = ORIG_DIR + "/out_0.0" + "_" + str(N) + "_" + str(SAMPLE_SIZE) + "_" + SUFFIX + ".het"
        dn = pd.read_csv(filen, sep = "\t", header = None,
                   names = ["DISPERSAL", "SELECTION", "WINDOW", "HET_AVG", "HET_STD"])
        dn_subset = dn[(dn["DISPERSAL"] == DISPERSAL) & (dn["WINDOW"] == 100000)]
        axs[sampling_idx,idx+3].axhline(dn_subset["HET_AVG"].iloc[0], linestyle = "--",color = "grey")#,label = "neutral")
        
        for color_idx,window in enumerate(window_sizes):
            LABEL = f"{window / 10**5:.1f}" + r"$\times 10^5$"
            d0_subset = d0[(d0["DISPERSAL"] == DISPERSAL) & (d0["WINDOW"] == window)]

            subset = d[(d["WINDOW"] == window)]
            subset = subset.sort_values(by='TIME')

            time =pd.concat((subset["TIME"], pd.Series(9999)))
            time = 10000 - time
            avg = pd.concat((subset["HET_AVG"], d0_subset["HET_AVG"]))
            std = pd.concat((subset["HET_STD"], d0_subset["HET_STD"]))
            axs[sampling_idx,idx+3].plot(time,avg, marker = "o", linestyle = "dotted", color = np.array(colors[sampling_idx])[color_idx],label = LABEL)
            axs[sampling_idx,idx+3].errorbar(time, avg, yerr=std, capsize=1, color = np.array(colors[sampling_idx])[color_idx], linestyle = "dotted")


        
pad = 2


LABELS = ["global sampling","local sampling"]
for ax, row in zip(axs[0:2,0], LABELS):
    ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                xycoords=ax.yaxis.label, textcoords='offset points', ha='right', va='center',rotation = 90)

for ax, row in zip(axs[0:2,-2], LABELS):
    ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                xycoords=ax.yaxis.label, textcoords='offset points', ha='right', va='center',rotation = 90)

col_headers = [r'$d$ = {:.3f}'.format(col) for col in [0.015,1, 0,0.015, 1]]
axes = fig.get_axes()
for ax in axes:
    sbs = ax.get_subplotspec()
    if sbs.is_first_row(): # column headers: neutral, s = 0.01, s = 0.1
        ax.annotate(col_headers[sbs.colspan.start], xy=(0.5, 1),xytext=(0, pad),
            xycoords="axes fraction",textcoords="offset points",ha="center",
            va="baseline")
        ax.get_xaxis().set_visible(False)
    if sbs.is_first_col():
        ax.set_ylabel("proportion of SNPs")
    else:
        ax.get_yaxis().set_visible(False)
    
    if sbs.colspan.start == 3:
        ax.set_ylabel("haplotype heterozygosity")
        ax.get_yaxis().set_visible(True)

    if(sbs.colspan.start < 2):
        ax.set_yscale('log')
        ax.set_ylim([1e-4,1e0])
        ax.set_xlabel("frequency")
        ax.set_xticks([0,0.5,1])
    else:
        ax.set_xscale("log")
        ax.set_ylim([0,1])
        ax.set_xticks(time)
        ax.set_xlabel("time since fixation")

handles, labels = axs[0,0].get_legend_handles_labels()
legend = fig.legend(handles, labels, bbox_to_anchor=(0.455, 0.88),markerscale=0.8, ncol = 2, columnspacing = 0.0005, handletextpad = 0.001, borderpad = 0.01, fontsize = 'small')#, title = "dispersal")

handles, labels = axs[0,3].get_legend_handles_labels()
legend = fig.legend(handles, labels, bbox_to_anchor=(0.815, 0.87),markerscale=0.8, ncol = 1, columnspacing = 0.001, handletextpad = 0.001, borderpad = 0.01, fontsize = 'small')#, title = "dispersal")

handles, labels = axs[1,3].get_legend_handles_labels()
legend = fig.legend(handles, labels, bbox_to_anchor=(0.815, 0.47),markerscale=0.8, ncol = 1, columnspacing = 0.001, handletextpad = 0.001, borderpad = 0.01, fontsize = 'small')#, title = "dispersal")

custom_line = Line2D([0], [0], color="grey", marker = "o", linestyle = "None", markersize=1)

# Add the legend with the custom line
fig.legend([custom_line],["neutral"], bbox_to_anchor=(0.45, 0.8), fontsize = 'small')

# fig.tight_layout()

for i, label in enumerate(('A', '', '','B', '')):
    axs[0,i].text(-0.35, 1.1, label, transform=axs[0,i].transAxes, fontsize = 8, weight = "bold")

# fig.tight_layout()
fig.savefig("figS9.pdf", bbox_inches = "tight", transparent = True)