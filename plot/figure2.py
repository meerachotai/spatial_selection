#!/usr/bin/env python3

import numpy as np
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

fig, ax = plt.subplots(nrows=1, ncols = 2, figsize = [3.2,2])

DIR = "data/"
SUBDIR = "tmrca"
disp_arr = [0.015, 0.02, 0.04, 0.1,0.5]

p = sns.diverging_palette(220, 20, s=150, l=45, n=len(disp_arr))
SAMPLE_SIZE = 100

OUT = "out_"
means = []
genome_length = 1e7
interval = 1e4
for idx,DISPERSAL_DISTANCE in enumerate(disp_arr):
    file = DIR + SUBDIR + "/" + OUT + str(DISPERSAL_DISTANCE) + ".tmrca"
    d = pd.read_csv(file, delim_whitespace = True, header = None,
                       names = ["TREE","LENGTH"] + [i for i in range(0,int(SAMPLE_SIZE/2))])

    bins = np.arange(0, genome_length, interval)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    d['window'] = pd.cut(np.cumsum(d['LENGTH']), bins = bins)
    print(d.shape)
    d.drop_duplicates(subset=['window'], keep = "first", inplace = True)
    print(d.shape)
    
    h_all = np.concatenate([np.tile(d.iloc[i,2:-1], (int(d['LENGTH'].iloc[i]), 1)) for i in range(len(d))]).flatten()
    if(DISPERSAL_DISTANCE == 0.5):
        DISPERSAL_DISTANCE = 1
    
    # chosen_trees = np.arange(0,len(d), 25)
    # h_all = np.concatenate([np.tile(d.iloc[i,2:], (int(d['LENGTH'][i]), 1)) for i in chosen_trees]).flatten()

    bins = np.quantile(h_all, [0, 0.25, 0.5, 0.75, 1]);  avg = np.mean(h_all)
    print(len(h_all),bins, avg)
    means.append(avg)
    sns.violinplot(y = h_all,x = DISPERSAL_DISTANCE, inner="quart", ax = ax[0],color = p[idx], linewidth = 0.2)
    
ax[0].set_ylabel("coalescence time")
ax[0].set_xlabel(r"dispersal ($d$)")

plt.setp(ax[0].xaxis.get_majorticklabels(), rotation=55, ha = "right",rotation_mode="anchor")
ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText = True)
ax[0].tick_params(axis='x', rotation=55)

# plt.setp(ax.collections, alpha=.3)
ax[0].scatter(x=range(len(means)),y=means,c="k")
disp_arr = [0.015, 0.02, 0.04, 0.1,1]
ax[0].set_xticklabels(["{:.3f}".format(DISPERSAL_DISTANCE) for DISPERSAL_DISTANCE in disp_arr])

# -----------------------------------Figure 2B-------------------------------------------------

DIR = "data/"
SUBDIR = "neutral_pi"

SUFFIXES = ["sampling_mid.", "sampling_rd.genome_", "sampling_rd."]
LABELS = ["local", "per genome", "global"]
PREFIX = "neutral"
N = 100 # number of replicates
SAMPLE_SIZE = 50
dispersal_array = [0.015, 0.02, 0.04, 0.1, 1.0]
cm = sns.diverging_palette(145, 300, s=60, n=len(SUFFIXES), center = "dark")

for idx,SUFFIX in enumerate(SUFFIXES):
    print(DIR + SUBDIR + "/" + PREFIX + "_" + SUFFIX + "pi")
    LABEL = LABELS[idx]
    if(LABEL == "per genome"):
        N = N * SAMPLE_SIZE
        print(DIR + SUBDIR + "/" + PREFIX + "_" + SUFFIX + "pi")
        d = pd.read_csv(DIR + SUBDIR + "/" + PREFIX + "_" + SUFFIX + "pi", delim_whitespace = True, header = None, names = ["DISPERSAL"] + [i for i in range(N)])
        d = d[d['DISPERSAL'] != 0.08]
        means = []
        N = 100
        for i in range(N):
            subd = d.iloc[:,(i*SAMPLE_SIZE)+1:((i+1)*SAMPLE_SIZE) + 1]
            means.append(subd.mean(axis = 1).tolist())
        avg = np.mean(means, axis = 0)
        std = np.std(means, axis = 0)
    else:
        N = 100
        d = pd.read_csv(DIR + SUBDIR + "/" + PREFIX + "_" + SUFFIX + "pi", delim_whitespace = True, header = None, names = ["DISPERSAL"] + [i for i in range(N)])
        d['DISPERSAL'] = d['DISPERSAL'].apply(lambda x: 1 if x == 1 else x)
        d = d[d['DISPERSAL'] != 0.08]
        d = d.sort_values(by='DISPERSAL')
        avg = d.iloc[:,1:].mean(axis = 1)
        std = d.iloc[:,1:].std(axis = 1)
    ax[1].plot(dispersal_array, avg, color = cm[idx], marker = "o", label = LABEL, linestyle = "dotted")
    ax[1].errorbar(dispersal_array, avg, yerr=std, fmt='o', capsize=2, color = cm[idx], linestyle = "dotted")

ax[1].set_xlabel(r"dispersal ($d$)")
ax[1].set_ylabel(r"heterozygosity ($\pi$)")
ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText = True)
ax[1].set_xscale("log")
ax[1].set_xticks(dispersal_array)

# ax[1].tick_params(axis='x', rotation=55)
ax[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.setp(ax[1].xaxis.get_majorticklabels(), rotation=60, ha = "right",rotation_mode="anchor") 

ax[1].set_ylim([0,7e-4])

expected = 4 * 10000 * 1e-8
ax[1].axhline(expected, linestyle = "--", color = "grey")

line1 = [Line2D([1], [0], color="grey", linestyle = "--")]
# labels1 = [r"$4N\mu$"]
# ax[1].add_artist(ax[1].legend(line1,labels1, bbox_to_anchor=(0.95, 0.25)))

handles, labels = ax[1].get_legend_handles_labels()
order = [2,0,1]
ax[1].legend([handles[idx] for idx in order],[labels[idx] for idx in order],bbox_to_anchor=(1, 0.4),loc='upper right')#, title = "sampling")


for i, label in enumerate(('a', 'b')):
    ax[i].text(-0.25, 1.1, label, transform=ax[i].transAxes, fontsize = 8,va='top', ha='right', weight = "bold")

    
fig.tight_layout()

fig.savefig("fig2.pdf", bbox_inches = "tight", transparent = True)