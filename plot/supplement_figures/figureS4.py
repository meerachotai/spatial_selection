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
matplotlib.rcParams['figure.dpi'] = 300

from matplotlib.lines import Line2D

gridspec = dict(wspace = 0.5)

DIR = "data/pfix_tfix/"
fig, ax = plt.subplots(nrows=1, ncols = 3, figsize = [4.6,1.5], gridspec_kw=gridspec)
generations = 500000
DOMINANCE = 0.5
N = 10000
cm = sns.color_palette(palette='Purples_d', n_colors = 3)


select_arr = [0.01,0.05, 0.1]
for idx, SELECTION in enumerate(select_arr):

    d = pd.read_csv(DIR + "/fixation_stats_" + str(SELECTION) + "_update.fix")
    
    d['dispersal'] = d['dispersal'].apply(lambda x: 1 if x == 0.5 else x)
    d = d[d['dispersal'] != 0.08]
    d["pfix_std"] = np.sqrt(2 * d["pfix"] * (1 - d["pfix"]) / d["total"])

    print(d["pfix"] * d["total"])
    expected_pfix = (1-np.exp(-2 * DOMINANCE * SELECTION))/(1-np.exp(-4 * N * DOMINANCE * SELECTION))
    
    y = np.array(d["pfix"])
    print(y)
    yerr = d["pfix_std"]

    ax[idx].plot(d["dispersal"], y, color = cm[idx], marker = "o", label = "{:.2f}".format(SELECTION), linestyle = "dotted")
    ax[idx].errorbar(d["dispersal"], y, yerr=yerr, fmt='o', capsize=2, color = np.array(cm)[idx])

    expected_pfix = (1-np.exp(-2 * DOMINANCE * SELECTION))/(1-np.exp(-4 * N * DOMINANCE * SELECTION))
    
    ax[idx].set_xlabel(r"dispersal ($d$)")
    ax[idx].axhline(expected_pfix, linestyle = "dashed", color = "grey")
    ax[idx].set_ylim([0,max(y) + (max(y) * 0.1)])

    ax[idx].set_title(f"$s = {SELECTION:.2f}$")
    
for a in ax:
    a.set_xscale("log")
    a.set_xticks(d["dispersal"])
    a.tick_params(axis='x', rotation=55)
    a.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.setp(a.xaxis.get_majorticklabels(), rotation=55, ha = "right",rotation_mode="anchor", fontsize = 4.8)

ax[0].set_ylabel("fixation probability")

fig.savefig("figS4.pdf",transparent=True, bbox_inches = "tight")#,pad_inches = 0)