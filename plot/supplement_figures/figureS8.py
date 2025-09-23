#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib
font = {'size': 5.5}
matplotlib.rc('font', **font)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['lines.linewidth'] = 0
matplotlib.rcParams['lines.markersize'] = 1
matplotlib.rcParams['figure.dpi'] = 300

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

fig, axs = plt.subplots(nrows=1, ncols = 2, figsize = [4,2])
bins = np.arange(0, 1.03, 0.04) 

# --------------------------------------Figure S8A----------------------------------------------------------------------------

disp_array = [0.05]

DIR = "data/supplement_data/garudh"
file = DIR + "/out_0.1_1000_100_500_sampling_rd_SGV_pop.garudh"
d = pd.read_csv(file, sep = "\t", header = None, index_col = False, names = ["DISPERSAL","SELECTION", "NSITES", "H1", "H12","H123", "H2/H1"])
colors = ["coral"]

for idx,DISPERSAL in enumerate(disp_array):
    h2h1 = d[d["DISPERSAL"] == DISPERSAL]["H2/H1"]
    hist, bins = np.histogram(h2h1, bins=bins)
    hist_cum = np.cumsum(hist) 
    axs[0].stairs(hist, bins, linestyle=(0, (1, 1)), color=colors[idx], label = "$x_0 = {:.5f}$".format(DISPERSAL))

disp_array = [0.015,0.5]

DIR = "data/supplement_data/garudh"
file = DIR + "/out_0.1_1000_100_500_sampling_rd.garudh"
d = pd.read_csv(file, sep = "\t", header = None, index_col = False, names = ["DISPERSAL","SELECTION", "NSITES", "H1", "H12","H123", "H2/H1"])
colors = ["indianred", "dodgerblue"]

for idx,DISPERSAL in enumerate(disp_array):
    h2h1 = d[d["DISPERSAL"] == DISPERSAL]["H2/H1"]
    hist, bins = np.histogram(h2h1, bins=bins)
    # bin_centers = 0.5 * (bins[1:] + bins[:-1])
    if(DISPERSAL == 0.5):
        DISPERSAL = 1
    axs[0].stairs(hist, bins, linestyle='-', color=colors[idx], label = "$d = {:.3f}$".format(DISPERSAL))
    
    # plt.hist(h2h1, label = "{:.3f}".format(DISPERSAL), bins = 50)

axs[0].set_xlabel("H2/H1")
# axs[0].set_ylim([0,250])
axs[0].set_xlim([0,1])
axs[0].legend(loc='upper left')

# --------------------------------------Figure S8B----------------------------------------------------------------------------

disp_array = [0.05]

DIR = "data/supplement_data/garudh"
file = DIR + "/out_0.1_1000_100_500_sampling_rd_SGV_pop.garudh"
d = pd.read_csv(file, sep = "\t", header = None, index_col = False, names = ["DISPERSAL","SELECTION", "NSITES", "H1", "H12","H123", "H2/H1"])
colors = ["coral"]

bins = np.arange(0.013,0.134, 0.004)

for idx,DISPERSAL in enumerate(disp_array):
    h2h1 = d[d["DISPERSAL"] == DISPERSAL]["H12"]
    hist, bins = np.histogram(h2h1, bins=bins)
    hist_cum = np.cumsum(hist) 
    # bin_centers = 0.5 * (bins[1:] + bins[:-1])
    axs[1].stairs(hist, bins, linestyle=(0, (1, 1)), color=colors[idx], label = "$x_0 = {:.5f}$".format(DISPERSAL))
    # plt.hist(h2h1, label = "{:.3f}".format(DISPERSAL), bins = 50)

# axs.set_xlabel("H12")
# axs.legend(loc='upper left')

disp_array = [0.015,0.5]

DIR = "data/supplement_data/garudh"
file = DIR + "/out_0.1_1000_100_500_sampling_rd.garudh"
d = pd.read_csv(file, sep = "\t", header = None, index_col = False, names = ["DISPERSAL","SELECTION", "NSITES", "H1", "H12","H123", "H2/H1"])
colors = ["indianred", "dodgerblue"]

for idx,DISPERSAL in enumerate(disp_array):
    h2h1 = d[d["DISPERSAL"] == DISPERSAL]["H12"]
    hist, bins = np.histogram(h2h1, bins=bins)
    # bin_centers = 0.5 * (bins[1:] + bins[:-1])
    if(DISPERSAL == 0.5):
        DISPERSAL = 1
    axs[1].stairs(hist, bins, linestyle='-', color=colors[idx], label = "$d = {:.3f}$".format(DISPERSAL))
    
    # plt.hist(h2h1, label = "{:.3f}".format(DISPERSAL), bins = 50)

axs[1].set_xlabel("H12")
# axs.set_ylim([0,350])
# axs.set_xlim([0,1])
# axs[1].legend(loc='upper right')

for i, label in enumerate(('a', 'b')):
    axs[i].text(-0.12, 1.05, label, transform=axs[i].transAxes, fontsize = 8, weight = "bold")

fig.tight_layout()
fig.savefig("figS8.pdf", bbox_inches = "tight", transparent = True)