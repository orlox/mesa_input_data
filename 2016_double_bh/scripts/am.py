#!/usr/bin/env python
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import rc
import mesa as ms
import math
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib.patheffects as path_effects
import os
import matplotlib.gridspec as gridspec
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size' : 20})
rc('text', usetex=True)
import numpy
from mesa import *
params = {'backend': 'pdf',
          'figure.figsize': [4.3, 5],
          'font.family':'serif',
          'font.size':10,
          'font.serif': 'Times Roman',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'medium',
          'legend.fontsize': 8,
          'legend.frameon' : False,
          'text.usetex': True,
          'figure.dpi': 600,
          'lines.markersize': 4,
          'lines.linewidth': 3,
          'lines.antialiased': False,
          'path.simplify': False,
          'legend.handlelength':3,
          'figure.subplot.bottom':0.1,
          'figure.subplot.top':0.975,
          'figure.subplot.left':0.15,
          'figure.subplot.right':0.95}

mpl.rcParams.update(params)

WHITE      = (1.00,1.00,1.00)
BLACK      = (0.00,0.00,0.00)
ORANGE     = (0.90,0.60,0.00)
SKY_BLUE   = (0.35,0.70,0.90)
BLUE_GREEN = (0.00,0.60,0.50)
YELLOW     = (0.95,0.90,0.25)
BLUE       = (0.00,0.45,0.70)
VERMILLION = (0.80,0.40,0.00)
RED_PURPLE = (0.80,0.60,0.70)

#hexcols[0] dark bluish
#hexcols[1] light blue
#hexcols[2] greenish
#hexcols[3] dark green
#hexcols[4] brownish
#hexcols[5] light brown
#hexcols[6] pinkish
#hexcols[7] dark something redish
#hexcols[8] magentish

hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',\
        '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466','#4477AA']

fig = plt.figure()

gs1 = gridspec.GridSpec(3, 1)
gs1.update(wspace =0, hspace = 0)
axes = []
axes.append(plt.subplot(gs1[0:1, :]))
axes.append(plt.subplot(gs1[1:2, :]))
axes.append(plt.subplot(gs1[2:3, :]))
folders = ["am_low", "am_mid", "am_high"]

mrange = np.arange(0.0001,100,0.1)
jsch = np.log10(4.6*1e16*mrange/3)
jkerr = np.log10(1.5*1e16*mrange/3)

for k in [0,1,2]:
    Z10 = history_data(folders[k], slname = "Z10.data", clean_starlog=False)
    Z20 = history_data(folders[k], slname = "Z20.data", clean_starlog=False)
    Z50 = history_data(folders[k], slname = "Z50.data", clean_starlog=False)

    axes[k].plot(Z50.get("mass"), Z50.get("log_j_rot"), label = "$Z=Z_\odot/50$", color = hexcols[1])
    axes[k].plot(Z20.get("mass"), Z20.get("log_j_rot"), label = "$Z=Z_\odot/20$", color = hexcols[5])
    axes[k].plot(Z10.get("mass"), Z10.get("log_j_rot"), label = "$Z=Z_\odot/10$", color = hexcols[3])

    axes[k].plot(mrange, jsch, label = "Schwarzchild", color = "k", ls = "--")
    axes[k].plot(mrange, jkerr, label = "Kerr", ls = ":", color = "k")

    if k==0:
        axes[k].legend(loc=4)

    axes[k].set_xlim(0,65)
    axes[k].set_ylim(15.6,18.4)

axes[0].text(0.08, 0.08, "$M_\mathrm{i}\simeq 50M_\odot,\;P_\mathrm{i}=0.9\;\mathrm{d}$", fontsize = 13, transform=axes[0].transAxes)
axes[1].text(0.08, 0.08, "$M_\mathrm{i}\simeq 63M_\odot,\;P_\mathrm{i}=1.0\;\mathrm{d}$", fontsize = 13, transform=axes[1].transAxes)
axes[2].text(0.08, 0.08, "$M_\mathrm{i}\simeq 80M_\odot,\;P_\mathrm{i}=1.1\;\mathrm{d}$", fontsize = 13, transform=axes[2].transAxes)

axes[2].set_xlabel("$m/M_\odot$")
axes[2].set_ylabel("$\log\;j_\mathrm{rot}\;\mathrm{[cm^2\;s^{-1}]}$")
axes[1].set_ylabel("$\log\;j_\mathrm{rot}\;\mathrm{[cm^2\;s^{-1}]}$")
axes[0].set_ylabel("$\log\;j_\mathrm{rot}\;\mathrm{[cm^2\;s^{-1}]}$")

axes[0].set_xticklabels([])
axes[1].set_xticklabels([])

plt.savefig("am.pdf", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                        transparent=False)#, bbox_inches='tight', pad_inches=0.0)

