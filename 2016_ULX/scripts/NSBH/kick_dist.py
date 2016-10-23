#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import *
import numpy as np
import sys
sys.path.insert(0, '../')
import kicks
from scipy.stats import maxwell
params = {'backend': 'pdf',
          'figure.figsize': [4.3, 2.2],
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
          'figure.subplot.bottom':0.2,
          'figure.subplot.top':0.95,
          'figure.subplot.left':0.15,
          'figure.subplot.right':0.92}

hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',\
        '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466','#4477AA']

mpl.rcParams.update(params)

A=np.array([np.append([vkick],kicks.sample_kick_distribution_P(23,5.5,55,1.4,vdist=lambda x:[float(vkick)], num_v=5, num_theta=400,num_phi=100)) for vkick in range(0,701,5)])

print(A)
print(A[:,0])
print(A[:,1])
print(A[:,2])

fig, axes= plt.subplots(1)

maxw = axes.fill_between(A[:,0],0,maxwell.pdf(A[:,0], scale=265.)/max(maxwell.pdf(A[:,0],scale=265.)),color="b", alpha=0.2, label="Maxwellian, $\\sigma=265~\\rm km~s^{-1}$")
merge, = axes.plot(A[:,0],10*A[:,1], color=hexcols[2],label="GW merge fraction $\\times$ 10")
disrupt, = axes.plot(A[:,0],A[:,2], color=hexcols[8],ls="--", label="Disrupt fraction")

axes.set_xlabel("$v_{\\rm kick}~\\rm[km~s^{-1}]$")
axes.set_ylabel("fraction")
#axes.set_xlim([0,50])
axes.set_ylim([0,1.19])
axes.legend([maxw,merge,disrupt],["Maxwellian, $\\sigma=265~\\rm km~s^{-1}$", "GW merge fraction $\\times$ 10", "Disrupt fraction"], loc="upper left", fontsize=7)

plt.savefig("kick_dist.pdf")
#plt.clf()
#plt.close(plt.gcf())
