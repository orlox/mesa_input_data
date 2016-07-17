#!/usr/bin/env python
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import rc
import mesa as ms
import math
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size' : 20})
rc('text', usetex=True)
import numpy
params = {'backend': 'pdf',
          'figure.figsize': [4.3, 3.5],
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
          'figure.subplot.bottom':0.15,
          'figure.subplot.top':0.95,
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

hist1 = ms.history_data("double_bh_Z50/1.9_0.8_1.10/LOGS1", slname = "history.data")
hist2 = ms.history_data("double_bh_Z50/1.9_0.8_1.10/LOGS2", slname = "history.data")

plt.plot(hist1.get("log_Teff")[279:],hist1.get("log_L")[279:], '-', color = hexcols[0])
plt.plot(hist2.get("log_Teff")[279:],hist2.get("log_L")[279:], '-', color = hexcols[5], lw=1.5)
plt.gca().text(0.12, 0.9, "$M_1=79M_\odot,M_2=64M_\odot,\;P_\mathrm{i}=1.1\;\mathrm{d}$", fontsize = 13, transform=plt.gca().transAxes)

plt.legend(loc=2)

plt.xlim([5.25,4.68])
plt.ylim([5.65,6.55])

plt.xlabel("$\log\;T_{\mathrm{eff}}$")
plt.ylabel("$\log\;L/L_\odot$")

ax = plt.gca()

axins = zoomed_inset_axes(ax, 2, loc=4,bbox_to_anchor=(0.6, 0.2), bbox_transform=ax.figure.transFigure)
plt.plot(hist1.get("log_Teff")[279:],hist1.get("log_L")[279:], '-', color = hexcols[0])
plt.plot(hist2.get("log_Teff")[279:],hist2.get("log_L")[279:], '-', color = hexcols[5], lw=1.5)
# sub region of the original image
x1, x2, y1, y2 = 4.78, 4.70, 5.682, 5.99
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticks([])
axins.set_yticks([])
axins.yaxis.tick_right()

mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")

fig = plt.gcf()

plt.gca().minorticks_on()

plt.savefig("HR_spiral.pdf", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                        transparent=False)#, bbox_inches='tight', pad_inches=0.0)

