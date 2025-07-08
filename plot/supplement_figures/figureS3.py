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

font = {'size': 4}
matplotlib.rc('font', **font)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['lines.linewidth'] = 0.7
matplotlib.rcParams['lines.markersize'] = 2
matplotlib.rcParams['figure.dpi'] = 300

from matplotlib.lines import Line2D

from matplotlib.gridspec import GridSpec

fig = plt.figure(figsize = [6,1.6])
gs = GridSpec(2, 7, width_ratios=[1, 0.05,1, 0.05,1, 0.005, 1])
gs.update(hspace=0.45)
ax1 = fig.add_subplot(gs[:, 0])
ax2 = fig.add_subplot(gs[:, 2])
ax3 = fig.add_subplot(gs[0, 4])
ax4 = fig.add_subplot(gs[0, 6])
ax5 = fig.add_subplot(gs[1, 4])
ax6 = fig.add_subplot(gs[1, 6])

pca_axes = [ax3,ax4,ax5,ax6]
# ----------------------------Figure S3A----------------------------
DIR = "data/supplement_data/fst_pca/"
file = DIR + "neutral.fst"

d = pd.read_csv(file, delim_whitespace = True, index_col = None, header = None)
d.rename(columns={0:'DISPERSAL'}, inplace=True)


d['DISPERSAL'] = d['DISPERSAL'].apply(lambda x: 1 if x == 0.5 else x)
d = d[d['DISPERSAL'] != 0.08]

d = d.sort_values(by='DISPERSAL')
avg = d.iloc[:,1:].mean(axis = 1)
std = d.iloc[:,1:].std(axis = 1)

dispersal_array = d['DISPERSAL']

ax2.errorbar(dispersal_array, avg, yerr=std, fmt='o', capsize=2, linestyle = "dotted")

ax2.set_xlabel(r"dispersal ($d$)")
ax2.set_ylabel(r"$F_{ST}$", labelpad = 0.5)

# ax2.ticklabel_format(axis='y', style='sci', scilimits=(-1,-1), useMathText = True)
ax2.set_xscale("log")
ax2.set_xticks(dispersal_array)
ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.setp(ax2.xaxis.get_majorticklabels(), rotation=60, ha = "right",rotation_mode="anchor") 

ax2.text(-0.1, 1.1, 'B', transform=ax2.transAxes, fontsize = 5,va='top', ha='right', weight = "bold")
# ----------------------------Figure S3B----------------------------

colors = ["indianred", "dodgerblue"]
DISPERSAL = 0.015
REPLICATE = 1
file = DIR + "/neutral_" + str(DISPERSAL) + "_" + str(REPLICATE) + ".loc"
d = pd.read_csv(file, sep = ",", index_col = 0, header = 0)
col = np.array(colors)[np.array(d['label'], dtype = int)]
ax1.scatter(d['x'], d['y'], marker='o', color = col, s = 0.3)
ax1.set_xlabel("x")
ax1.set_ylim(0,1)
ax1.set_xlim(0,1)
ax1.set_ylabel("y")

ax1.text(-0.1, 1.1, 'A', transform=ax1.transAxes, fontsize = 5,va='top', ha='right', weight = "bold")
# ----------------------------Figure S3C----------------------------

dispersal_array = [0.015, 0.04, 0.1, 0.5]
REPLICATE = 1
colors = ["indianred", "dodgerblue"]
for DISPERSAL,axs in zip(dispersal_array, pca_axes):
    
    file = DIR + "/neutral_" + str(DISPERSAL) + "_" + str(REPLICATE) + ".pca"
    d = pd.read_csv(file, sep = ",", index_col = 0, header = 0)
    col = np.array(colors)[np.array(d['label'], dtype = int)]
    axs.scatter(d['PC1'], d['PC2'], marker='o', color = col, s = 0.3)
    if(DISPERSAL == 0.5):
        DISPERSAL = 1.0
    axs.set_title(r"$d = {:.3f}$".format(DISPERSAL), pad = 1)


pca_axes[0].set_ylabel("PC2", labelpad=0.4)
pca_axes[2].set_ylabel("PC2", labelpad=0.4)

pca_axes[2].set_xlabel("PC1")
pca_axes[3].set_xlabel("PC1")
pca_axes[0].text(-0.1, 1.22, 'C', transform=pca_axes[0].transAxes, fontsize = 5,va='top', ha='right', weight = "bold")

fig.tight_layout()
fig.savefig("figS3.pdf",bbox_inches = "tight", transparent = True)