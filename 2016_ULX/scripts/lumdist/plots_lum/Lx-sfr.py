#!/usr/bin/env python
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import matplotlib.patheffects as pe
import matplotlib.gridspec as grd
import matplotlib as mpl
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
          'lines.markersize': 2,
          'lines.linewidth': 3,
          'lines.antialiased': False,
          'path.simplify': False,
          'legend.handlelength':3,
          'figure.subplot.bottom':0.08,
          'figure.subplot.top':0.95,
          'figure.subplot.left':0.15,
          'figure.subplot.right':0.9}

hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',\
        '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466','#4477AA']

mpl.rcParams.update(params)

Lsun = 3.9e33

plt.figure()
gs = grd.GridSpec(2,2, width_ratios=[25,1], wspace=0.05, hspace=0.1)
axes0 = plt.subplot(gs[2])
axes1 = plt.subplot(gs[0])

axescb0 = plt.subplot(gs[3])
axescb1 = plt.subplot(gs[1])

inputs = ["lumdist_-2.500.dat","lumdist_-4.000.dat"]
titles = ["$\\log Z=-2.5$","$\\log Z=-4$"]
for i in [0,1]:
    data = np.loadtxt(inputs[i], skiprows=1)
    if i == 1:
        axes = axes0
    else:
        axes = axes1

    data = np.clip(data,1e-10,1e10)

    axes.fill_between(np.log10(data[:,0]),np.log10(data[:,5]),np.log10(data[:,6]),\
            color=mpl.colors.colorConverter.to_rgba(hexcols[0], alpha = 0.1), lw=0)
    axes.fill_between(np.log10(data[:,0]),np.log10(data[:,6]),np.log10(data[:,7]),\
            color=mpl.colors.colorConverter.to_rgba(hexcols[0], alpha = 0.35), lw=0)
    axes.fill_between(np.log10(data[:,0]),np.log10(data[:,7]),np.log10(data[:,8]),\
            color=mpl.colors.colorConverter.to_rgba(hexcols[0], alpha = 0.6), lw=0)
    axes.fill_between(np.log10(data[:,0]),np.log10(data[:,8]),np.log10(data[:,10]),\
            color=mpl.colors.colorConverter.to_rgba(hexcols[0], alpha = 0.85), lw=0)
    axes.fill_between(np.log10(data[:,0]),np.log10(data[:,10]),np.log10(data[:,11]),\
            color=mpl.colors.colorConverter.to_rgba(hexcols[0], alpha = 0.6), lw=0)
    axes.fill_between(np.log10(data[:,0]),np.log10(data[:,11]),np.log10(data[:,12]),\
            color=mpl.colors.colorConverter.to_rgba(hexcols[0], alpha = 0.35), lw=0)
    axes.fill_between(np.log10(data[:,0]),np.log10(data[:,12]),np.log10(data[:,13]),\
            color=mpl.colors.colorConverter.to_rgba(hexcols[0], alpha = 0.1), lw=0)

    sfra = np.linspace(0.1,4,50)
    valsa = 2.6*np.power(sfra,1.0/0.6)
    sfrb = np.linspace(4,100,50)
    valsb = 6.7*sfrb
    axes.plot(np.log10(np.append(sfra,sfrb)), np.log10(np.append(valsa,valsb)), hexcols[6], label="Grimm et al. (2003)")

    axes.plot(np.log10(data[:,0]), np.log10(data[:,9]), hexcols[5], label="Median",ls="-")
    axes.plot(np.log10(data[:,0]), np.log10(data[:,1]), hexcols[2], label="Linear",ls="--")

    axes.legend(loc="upper left", title = titles[i])

#axes1.text(1, 1.0, '$1^{\\rm st}$ decile',
#         horizontalalignment='center',
#         verticalalignment='top',
#         multialignment='center', fontsize=10, color=hexcols[0], alpha=0.5, rotation = 35)
#axes1.text(-0.2, 1.3, '$9^{\\rm th}$ decile',
#         horizontalalignment='center',
#         verticalalignment='top',
#         multialignment='center', fontsize=10, color=hexcols[0], alpha=0.5, rotation = 25)

axes1.set_xlim([-1,2])
axes1.set_ylim([0,3.8])
axes0.set_xlim([-1,2])
axes0.set_ylim([0,3.8])
axes1.xaxis.set_ticklabels([])

axes0.set_xlabel("$\\log {\\rm SFR}\;{\\rm [M_\\odot\\;yr^{-1}]}$")
axes1.set_ylabel("$\\log L_{\\rm X}\;{\\rm[10^{39}\\; erg\;s^{-1}]}$")
axes0.set_ylabel("$\\log L_{\\rm X}\;{\\rm[10^{39}\\; erg\;s^{-1}]}$")

cmap = mpl.colors.ListedColormap([mpl.colors.colorConverter.to_rgba(hexcols[0], alpha = 0.85),\
                                  mpl.colors.colorConverter.to_rgba(hexcols[0], alpha = 0.6),\
                                  mpl.colors.colorConverter.to_rgba(hexcols[0], alpha = 0.35),\
                                  mpl.colors.colorConverter.to_rgba(hexcols[0], alpha = 0.1),\
                                  "w"])
bounds = [0, 20, 40, 60, 80, 100]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
cbar = mpl.colorbar.ColorbarBase(axescb0, cmap=cmap,
                                norm=norm,
                                ticks=bounds,
                                spacing='proportional',
                                orientation='vertical')
cbar.ax.set_yticklabels(['$0\\%$','$20\\%$','$40\\%$','$60\\%$','$80\\%$','$100\\%$'])
cbar = mpl.colorbar.ColorbarBase(axescb1, cmap=cmap,
                                norm=norm,
                                ticks=bounds,
                                spacing='proportional',
                                orientation='vertical')
cbar.ax.set_yticklabels(['$0\\%$','$20\\%$','$40\\%$','$60\\%$','$80\\%$','$100\\%$'])

plt.savefig("Lx-sfr.pdf")
