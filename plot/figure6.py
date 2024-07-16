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

gridspec = dict(hspace=0.1, height_ratios=[1, 1, 0.3, 1], wspace = 0.05)

fig, axs = plt.subplots(nrows=4, ncols = 3, figsize = [3.6,5], gridspec_kw=gridspec)
for i in range(3):
    axs[2,i].set_visible(False)

# ----------------------------Figure 6A---------------------------------------------------------------

DIR = "data/"
SUBDIR = "sfs"

PREFIX = "sfs_window"
SAMPLE_SIZE = 100
SUFFIXES = ["sampling_rd","sampling_mid"]
LABELS = ["global sampling","local sampling"]

colors = sns.diverging_palette(220, 20, s=150, l=45, n=6)
x = np.arange(0, SAMPLE_SIZE + 1, 1)
disp_arr = [0.015, 0.02, 0.04, 0.1,0.5]

select_arr = [0,0.01, 0.1]
window_sizes = [1e7,150000,1000000]
window_sizes = [float(win) for win in window_sizes]
select_arr = [float(s) for s in select_arr]

for sampling_idx, SUFFIX in enumerate(SUFFIXES):
    for idx,SELECTION in enumerate(select_arr):
        window_size = window_sizes[idx]
        d = pd.read_csv(DIR + SUBDIR + "/" + PREFIX + "_" + str(SELECTION) + "_" + str(window_size) + "_" + SUFFIX + ".sfs", delim_whitespace = True, header = None,
                   names = ["TYPE","DISPERSAL", "SELECTION"] + [i for i in range(SAMPLE_SIZE + 1)] + ["FILE_COUNT"])
        if(SELECTION == 0):
            print(SELECTION, SUFFIX)
#             print(d.iloc[:,:10][d["TYPE"] == "AVG"])
            print(d.iloc[:,45:55][d["TYPE"] == "AVG"])
        for color_idx,DISPERSAL_DISTANCE in enumerate(disp_arr):
            
            sfs = d[(d["DISPERSAL"] == DISPERSAL_DISTANCE) & (d["TYPE"] == "AVG") & (d["SELECTION"] == SELECTION)].iloc[:, 3:]
            if(DISPERSAL_DISTANCE == 0.5):
                DISPERSAL_DISTANCE = 1
            
            LABEL = r"$d$ = " + '{:.3f}'.format(DISPERSAL_DISTANCE)
            
            if(sampling_idx == 0):
                axs[sampling_idx,idx].scatter(sfs.columns[1:-1]/100,sfs.iloc[0,1:-1].tolist(), color = np.array(colors)[color_idx], label = LABEL, s = 0.8,alpha = 0.8)
            else:
                axs[sampling_idx,idx].scatter(sfs.columns[1:-1]/100,sfs.iloc[0,1:-1].tolist(), color = np.array(colors)[color_idx], s = 0.8, alpha = 0.8)

# ---------------------------Figure GB----------------------------------------------------------------

LABELS = ["global","local"]
# cm = sns.diverging_palette(220, 20, s=150, l=45, n=len(SUFFIXES))
cm = sns.diverging_palette(300,145, s=60, n = len(SUFFIXES))

for sampling_idx, SUFFIX in enumerate(SUFFIXES):
    for idx,SELECTION in enumerate(select_arr):
        window_size = window_sizes[idx]
        d = pd.read_csv(DIR + SUBDIR + "/" + PREFIX + "_" + str(SELECTION) + "_" + str(window_size) + "_" + SUFFIX + ".tajimad", delim_whitespace = True, header = None,
               names = ["DISPERSAL", "SELECTION", "WINDOW","AVG", "STD", "NFILES"])
        d['DISPERSAL'] = d['DISPERSAL'].apply(lambda x: 1 if x == 0.5 else x)
        d = d[d['DISPERSAL'] != 0.08]
#         print(d)
        d_sorted = d.sort_values(by='DISPERSAL') # sorting by dispersal to add dispersal = 1 to end
        axs[3, idx].errorbar(d_sorted["DISPERSAL"], d_sorted["AVG"], yerr = d_sorted["STD"], marker = "o", color = cm[sampling_idx], label = LABELS[sampling_idx],linestyle = "dotted", capsize=1)
        axs[3, idx].axhline(0, linestyle = "dashed", color = "grey")

# (_, caps, _) = plt.errorbar(
    # x_values, y_values, yerr=y_error, fmt='o', markersize=8, capsize=20)
# for cap in caps:
    # cap.set_markeredgewidth(1)
# -------------------------------------------------------------------------------------------

pad = 2
       
col_headers = [r'$s$ = {:.2f}'.format(col) if col != 0 else "neutral" for col in select_arr]

axes = fig.get_axes()
for ax in axes:
    sbs = ax.get_subplotspec()
    if sbs.is_first_row(): # column headers: neutral, s = 0.01, s = 0.1
        ax.annotate(col_headers[sbs.colspan.start], xy=(0.5, 1),xytext=(0, pad),
            xycoords="axes fraction",textcoords="offset points",ha="center",
            va="baseline")
        ax.get_xaxis().set_visible(False)
    elif sbs.rowspan.start == 3:
        ax.annotate(col_headers[sbs.colspan.start], xy=(0.5, 1),xytext=(0, pad),
            xycoords="axes fraction",textcoords="offset points",ha="center",
            va="baseline")
        
    # y-axis labels only on the first column
    if sbs.is_first_col():
        if sbs.rowspan.start < 2: # second row
            ax.set_ylabel("proportion of SNPs")
        else:
            ax.set_ylabel("Tajima's D")
    else:
        ax.get_yaxis().set_visible(False)
        
        
    # x-axis labels only on the bottom row
    # set x,y-axis logscale and limits
    
    if sbs.rowspan.start < 2: # second row
        ax.set_ylim([1e-4,1e0])
        ax.set_yscale('log')
        if sbs.rowspan.start == 1:
            ax.set_xticks([0,0.5,1])
            ax.set_xlabel("frequency")
    else:
        ax.set_ylim([-2,1])
        ax.set_xscale("log")
        ax.set_xticks(d_sorted["DISPERSAL"])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=55, ha = "right",rotation_mode="anchor", fontsize = 4.8)
        ax.set_yticks(np.arange(-2,3,1))
        ax.set_xlabel(r"dispersal ($d$)")

LABELS = ["global sampling","local sampling"]
for ax, row in zip(axs[0:2,0], LABELS):
    ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                xycoords=ax.yaxis.label, textcoords='offset points', ha='right', va='center',rotation = 90)
    
for i, label in enumerate(('A', '', '','B')):
    axs[i,0].text(-0.68, 1.1, label, transform=axs[i,0].transAxes, fontsize = 8, weight = "bold")

# ------------------------------------------------------------------

x = np.arange(0.01,1,0.01)
y = 1/x
summed_y = sum(y)
y = y / summed_y
axs[0,0].plot(x,y, linestyle = "--", color = "black")
axs[1,0].plot(x,y, linestyle = "--", label = "$1/x$", color = "black")


handles, labels = axs[0,0].get_legend_handles_labels()
legend = fig.legend(handles, labels, bbox_to_anchor=(0.38, 0.88),markerscale=0.8, ncol = 2, columnspacing = 0.001, handletextpad = 0.001, borderpad = 0.01, fontsize = 'small')#, title = "dispersal")

handles, labels = axs[1,0].get_legend_handles_labels()
legend = fig.legend(handles, labels, bbox_to_anchor=(0.38, 0.82),markerscale=1, columnspacing = 0.01, handletextpad = 0.01)

handles, labels = axs[3,0].get_legend_handles_labels()
legend = fig.legend(handles, labels, bbox_to_anchor=(0.37, 0.22), title = "sampling")


fig.savefig("fig6lefttop.pdf", bbox_inches = "tight", transparent = True)

# -------------------------Figure 6C-----------------------------------------

import math

DIR = "data/"
SUBDIR = "haplotype_frequency_spectra"
tothap = 100
SAMPLE_SIZE = 100
SUFFIXES = ["sampling_rd", "sampling_mid"]
LABEL = ["global", "local"]
SELECTION = 0.1

OUT = DIR + SUBDIR + "/out"
disp_arr = [0.5,0.015]

palette = sns.husl_palette(5, h = 0.5, s = 0.8)
palette = palette * math.ceil(tothap/len(palette))

fig, axs = plt.subplots(nrows=3, ncols = 1, figsize = [3.6,1], sharex = True)

plot_idx = 0

for idx, DISPERSAL_DISTANCE in enumerate(disp_arr):
    for sampling_idx, SUFFIX in enumerate(SUFFIXES):
        
        d = pd.read_csv(OUT + "_" + str(SELECTION) + "_" + str(DISPERSAL_DISTANCE) + "_" + SUFFIX + ".hfs", delim_whitespace = True, header = None, index_col = None,
                       names = ["DISPERSAL","SELECTION"] + [i for i in range(tothap)])
        
        meanhc = d.iloc[:,2:].mean(axis = 0) / SAMPLE_SIZE
        
        sns.despine(ax=axs[plot_idx], left=True)
    
        x1 = 0
        for i, c in enumerate(meanhc):
            x2 = x1 + c
            if c > (1/SAMPLE_SIZE):
                color = palette[i]
            else:
                color = "lightgrey"
            axs[plot_idx].axvspan(x1, x2, facecolor=color, edgecolor = "white")
            x1 = x2
        axs[plot_idx].set_yticks([])
        
        if(DISPERSAL_DISTANCE == 0.5):
            # label = "panmictic"
            label = "dispersal = " + "{:.3f}".format(DISPERSAL_DISTANCE * 2) + "\n" + LABEL[sampling_idx] + " sampling"
        else:
            label = "dispersal = " + "{:.3f}".format(DISPERSAL_DISTANCE) + "\n" + LABEL[sampling_idx] + " sampling"
        axs[plot_idx].set_ylabel(label, rotation=0, labelpad=40, loc = "bottom")
        plot_idx = plot_idx + 1
        if(idx == 0):
            break
for i, label in enumerate(('C', '')):
    axs[i].text(-0.21, 1.2, label, transform=axs[i].transAxes, fontsize = 8,va='top', ha='right', weight = "bold")

fig.savefig("fig6leftbottom.pdf", bbox_inches = "tight", transparent = True)

# ---------------------------Figure 6D---------------------------------------------------------

gridspec = dict(hspace=0.1, height_ratios=[1, 1, 0.3, 1, 1], wspace = 0.05)

fig, axs = plt.subplots(nrows=5, ncols = 3, figsize = [3.6,6.5], gridspec_kw=gridspec)
for i in range(3):
    axs[2,i].set_visible(False)

select_arr = [0.01,0.05, 0.1]

DOMINANCE = 0.5
TOTAL_FILES = 10
SAMPLE_SIZE = 100
N = 1000

MU=1e-8
SEED = 5
genome_length = 1e7
n = 2

DIR = "data/"
SUBDIR = "haplotype_statistics"
SUFFIXES = ["sampling_rd","sampling_mid"]
PREFIX = "out"

window_sizes = [10000, 30000,100000]

colors_all = sns.diverging_palette(145, 300, s=60, l = 50, n = len(SUFFIXES) * len(window_sizes))# + 2)
colors_local = colors_all[:3][::-1]
colors_global = colors_all[-3:]
colors = [colors_global] + [colors_local]

LABELS = ["global sampling","local sampling"]
import seaborn as sns

x = np.arange(0, SAMPLE_SIZE + 1, 1)
disp_arr = [0.015, 0.02, 0.04, 0.1, 0.5]
select_arr = [0,0.01,0.1]
select_arr = [float(s) for s in select_arr]

for sampling_idx, SUFFIX in enumerate(SUFFIXES):
    for idx,SELECTION in enumerate(select_arr):
        file = DIR + SUBDIR + "/" + PREFIX + "_" + str(SELECTION) + "_" + str(N) + "_" + str(SAMPLE_SIZE) + "_" + SUFFIX + ".het"
        
        d = pd.read_csv(file, delim_whitespace = True, header = None,
                   names = ["DISPERSAL", "SELECTION", "WINDOW", "HET_AVG", "HET_STD"])
        d['DISPERSAL'] = d['DISPERSAL'].apply(lambda x: 1 if x == 0.5 else x)
        d = d[d['DISPERSAL'] != 0.08]

        for color_idx,window in enumerate(window_sizes):
            # if(SELECTION == 0 and color_idx == 0):
            #     continue
            LABEL = f"{window / 10**5:.1f}" + r"$\times 10^5$"
            subset = d[(d["SELECTION"] == SELECTION) & (d["WINDOW"] == window)]
            subset = subset.sort_values(by='DISPERSAL')
            axs[sampling_idx,idx].plot(subset["DISPERSAL"],subset["HET_AVG"], marker = "o", linestyle = "dotted", color = np.array(colors[sampling_idx])[color_idx],label = LABEL)
            axs[sampling_idx,idx].errorbar(subset["DISPERSAL"],subset["HET_AVG"], yerr=subset["HET_STD"], capsize=1, color = np.array(colors[sampling_idx])[color_idx], linestyle = "dotted")

    
# handles, labels = axs[1,1].get_legend_handles_labels()
# legend = fig.legend(handles, labels, bbox_to_anchor=(1.15, 0.65),markerscale=1.25, title = "window size")


# ---------------------------Figure 6E---------------------------------------------------------

window_sizes = [10000, 30000,100000]
for sampling_idx, SUFFIX in enumerate(SUFFIXES):
    for idx,SELECTION in enumerate(select_arr):
        file = DIR + SUBDIR + "/" + PREFIX + "_" + str(SELECTION) + "_" + str(N) + "_" + str(SAMPLE_SIZE) + "_" + SUFFIX + ".nhaplo"
    
        d = pd.read_csv(file, delim_whitespace = True, header = None,
                   names = ["DISPERSAL", "SELECTION", "WINDOW", "HET_AVG", "HET_STD"])
        d['DISPERSAL'] = d['DISPERSAL'].apply(lambda x: 1 if x == 0.5 else x)
        d = d[d['DISPERSAL'] != 0.08]
        
        for color_idx,window in enumerate(window_sizes):
            # if(SELECTION == 0 and color_idx == 0):
            #     continue
            LABEL = f"{window / 10**5:.1f}" + r"$\times 10^5$"
            subset = d[(d["SELECTION"] == SELECTION) & (d["WINDOW"] == window)]
            subset = subset.sort_values(by='DISPERSAL')
            axs[sampling_idx + 3,idx].plot(subset["DISPERSAL"],subset["HET_AVG"], marker = "o", linestyle = "dotted", color = np.array(colors[sampling_idx])[color_idx],label = LABEL)
            axs[sampling_idx + 3,idx].errorbar(subset["DISPERSAL"],subset["HET_AVG"], yerr=subset["HET_STD"], capsize=1, color = np.array(colors[sampling_idx])[color_idx],linestyle = "dotted")


pad = 2

col_headers = [r'$s$ = {:.2f}'.format(col) if col != 0 else "neutral" for col in select_arr]
axes = fig.get_axes()
for ax in axes:
    sbs = ax.get_subplotspec()
    if sbs.is_first_row() or sbs.rowspan.start == 3:
        ax.annotate(col_headers[sbs.colspan.start], xy=(0.5, 1),xytext=(0, pad),
            xycoords="axes fraction",textcoords="offset points",ha="center",
            va="baseline")
    
    if sbs.rowspan.start == 1 or sbs.rowspan.start == 4:
        ax.set_xlabel(r"dispersal ($d$)")
    else:
        ax.get_xaxis().set_visible(False)
        
    if(sbs.rowspan.start < 2):
        ax.set_ylim([0,1])
    else:
        ax.set_ylim([0,80])
        
    if sbs.is_first_col():
        if(sbs.rowspan.start < 2):
            ax.set_ylabel("haplotype heterozygosity") 
        else:
            ax.set_ylabel("number of haplotypes")
    else:
        ax.get_yaxis().set_visible(False)
    
for ax, row in zip(axs[0:2,0], LABELS):
    ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                xycoords=ax.yaxis.label, textcoords='offset points',
                ha='right', va='center',rotation = 90)

    
for ax, row in zip(axs[3:,0], LABELS):
    ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                xycoords=ax.yaxis.label, textcoords='offset points',
                ha='right', va='center',rotation = 90)

axes = axs.flatten()
for ax in axes:
    ax.set_xscale("log")
    ax.set_xticks(subset["DISPERSAL"])
    ax.tick_params(axis='x', rotation=55)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=55, ha = "right",rotation_mode="anchor", fontsize = 4.6)

handles, labels = axs[0,1].get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor=(0.28, 0.80),markerscale=1, handletextpad = 0.01, title = "window size", labelspacing = 0.2, borderpad = 0.05)

handles, labels = axs[1,1].get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor=(0.28, 0.62),markerscale=1, handletextpad = 0.01, title = "window size", labelspacing = 0.2, borderpad = 0.05)

handles, labels = axs[-2,1].get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor=(0.88, 0.46),markerscale=1, handletextpad = 0.01, title = "window size", labelspacing = 0.2, borderpad = 0.05)

handles, labels = axs[-1,1].get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor=(0.88, 0.27),markerscale=1, handletextpad = 0.01, title = "window size", labelspacing = 0.2, borderpad = 0.05)

for i, label in enumerate(('D', '', '','E', '')):
    axs[i,0].text(-0.65, 1.1, label, transform=axs[i,0].transAxes, fontsize = 8, weight = "bold")

fig.savefig("fig6right.pdf", bbox_inches = "tight", transparent = True)




