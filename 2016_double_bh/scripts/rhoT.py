#!/usr/bin/env python
from mesa import *
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
import numpy
import matplotlib.patches as patches
from scipy.interpolate import griddata
#from matplotlib.mlab import griddata
import math
import scipy
from matplotlib import ticker

params = {'backend': 'pdf',
          'figure.figsize': [4.3, 3.0],
          'font.family':'serif',
          'font.size':10,
          'font.serif': 'Times Roman',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'medium',
          'legend.fontsize': 8,
          'legend.frameon' : False,
          'text.usetex': True,
          'figure.dpi': 600,
          'lines.markersize': 2,
          'lines.linewidth': 3,
          'lines.antialiased': False,
          'path.simplify': False,
          'legend.handlelength':3,
          'figure.subplot.bottom':0.15,
          'figure.subplot.top':0.9,
          'figure.subplot.left':0.15,
          'figure.subplot.right':0.9}

mpl.rcParams.update(params)

hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',\
        '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466','#4477AA']


pur = numpy.loadtxt("pur.data", usecols=[0,1])

high = history_data('.', slname = "high.data", clean_starlog=False)
pisn = history_data('.', slname = "pisn.data", clean_starlog=False)
lowm = history_data('.', slname = "lowm.data", clean_starlog=False)

plt.plot(high.get("log_center_Rho"), high.get("log_center_T"), label = "$\sim 200 M_\odot$", color = hexcols[1])
plt.plot(pisn.get("log_center_Rho"), pisn.get("log_center_T"), label = "$\sim 90 M_\odot$", color = hexcols[3])
plt.plot(lowm.get("log_center_Rho"), lowm.get("log_center_T"), label = "$\sim 35 M_\odot$", color = hexcols[5])

plt.legend(loc=4)

plt.xlim([2.5,8.1])
plt.ylim([8.5,9.85])

vert_x = pur[:,0]
vert_y = pur[:,1]

p = plt.Polygon(np.column_stack((vert_x, vert_y)), facecolor='gray', alpha=.45, edgecolor='none')

text(2.8, 9.1, "$e^+e^-$\\\\\\\\$\Gamma<4/3$", fontsize = 15, color = "k")

plt.gca().add_artist(p)

plt.ylabel("$\log\;T_{\mathrm{center}}\;\mathrm{[K]}$")
plt.xlabel("$\log\;\\rho_\mathrm{center}\;\mathrm{[g\;cm^{-3}]}$")

plt.savefig("rhoT.pdf", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                        transparent=False)
