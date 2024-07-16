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
matplotlib.rcParams['lines.markersize'] = 1

from matplotlib.lines import Line2D

gridspec = dict(wspace = 0.5)

DIR = "data/"
SUBDIR = "pfix_tfix"
fig, ax = plt.subplots(nrows=1, ncols = 2, figsize = [3.2,1.5], gridspec_kw=gridspec)
generations = 500000
DOMINANCE = 0.5
N = 10000
cm = sns.color_palette(palette='Purples_d', n_colors = 3)


select_arr = [0.01,0.05, 0.1]
for idx, SELECTION in enumerate(select_arr):
    print(DIR + "/" + SUBDIR + "/fixation_stats_" + str(SELECTION) + ".fix")
    d = pd.read_csv(DIR + "/" + SUBDIR + "/fixation_stats_" + str(SELECTION) + ".fix")
    d['dispersal'] = d['dispersal'].apply(lambda x: 1 if x == 0.5 else x)
    d = d[d['dispersal'] != 0.08]
    d["pfix_std"] = np.sqrt(d["pfix"] * (1 - d["pfix"]) / d["total"])
    print(d)

    print(d["pfix"] * d["total"])
    expected_pfix = (1-np.exp(-2 * DOMINANCE * SELECTION))/(1-np.exp(-4 * N * DOMINANCE * SELECTION))
    # pfix_panmictic = float(d[d["dispersal"] == 1]["pfix"])
    
    y = np.array(d["pfix"]/ expected_pfix)
    print(y)
    yerr = d["pfix_std"] / expected_pfix

    ax[1].plot(d["dispersal"], y, color = cm[idx], marker = "o", label = "{:.2f}".format(SELECTION), linestyle = "dotted")
    ax[1].errorbar(d["dispersal"], y, yerr=yerr, fmt='o', capsize=2, color = np.array(cm)[idx])

    tfix_panmictic = float(d[d["dispersal"] == 1]["tfix_mean"]) - generations
    
    y = np.array((d["tfix_mean"] - generations)/ tfix_panmictic)
    print(y)
    yerr = d["tfix_std"] / tfix_panmictic

    ax[0].plot(d["dispersal"][:-1], y[:-1], color = cm[idx], marker = "o", label = "s = {:.2f}".format(SELECTION), linestyle = "dotted")
    ax[0].errorbar(d["dispersal"][:-1], y[:-1], yerr=yerr[:-1], fmt='o', capsize=2, color = np.array(cm)[idx]) # remove the d = 1

# matplotlib.rcParams['xtick.minor.size'] = 0
# matplotlib.rcParams['xtick.minor.width'] = 0


ax[0].set_xscale("log")
ax[0].xaxis.set_major_locator(matplotlib.ticker.NullLocator())
ax[0].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
print(d["dispersal"][:-1])
ax[0].set_xticks(d["dispersal"][:-1], labels = d["dispersal"][:-1])
print("Ticks before:",list(ax[0].get_xticklabels()))

ax[0].tick_params(axis='x', rotation=55)
ax[0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.setp(ax[0].xaxis.get_majorticklabels(), rotation=55, ha = "right",rotation_mode="anchor", fontsize = 4.8) 

ax[1].set_xscale("log")
ax[1].set_xticks(d["dispersal"])
ax[1].tick_params(axis='x', rotation=55)
ax[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.setp(ax[1].xaxis.get_majorticklabels(), rotation=55, ha = "right",rotation_mode="anchor", fontsize = 4.8)

for a in ax:
    a.set_xlabel(r"dispersal ($d$)")
    a.axhline(1, linestyle = "dashed", color = "grey")

ax[1].set_yticks(np.arange(0.5,1.3,0.25))
ax[0].set_yticks(np.arange(0,2.6,0.5))

ax[1].set_ylabel("fixation probability (relative to expectation)")
ax[0].set_ylabel("time to fixation (relative to panmictic)")

lines,labels = ax[1].get_legend_handles_labels()
ax[1].legend(lines, labels,bbox_to_anchor=(0.75, 0.45), loc='upper right',markerscale=0.8)

for i, label in enumerate(('A', 'B')):
    ax[i].text(-0.4, 1.2, label, transform=ax[i].transAxes, fontsize = 8,va='top', ha='right', weight = "bold")
    
# fig.tight_layout()
fig.savefig("fig4.pdf",transparent=True, bbox_inches = "tight")#,pad_inches = 0)
