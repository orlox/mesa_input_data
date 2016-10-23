#!/usr/bin/env python
from mesa_data import *
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import matplotlib.patheffects as pe
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
          'figure.subplot.right':0.92}

hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',\
        '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466','#4477AA']

mpl.rcParams.update(params)

fig, axes= plt.subplots(1)

profs = [Mesa_Data("binary_history_250.data"), Mesa_Data("binary_history_300.data"), Mesa_Data("binary_history_350.data")]
colors = [hexcols[2], hexcols[8], hexcols[5]]
labels = ["$\\log Z=-2.5, P_{\\rm orb,i}=1.4{\\rm d}$",\
        "$\\log Z=-3.0, P_{\\rm orb,i}=1.2{\\rm d}$",\
        "$\\log Z=-3.5, P_{\\rm orb,i}=1.1{\\rm d}$"]
for i in [0,1,2]:
    axes.plot(profs[i].get("age")/1e6,np.log10(profs[i].get("period_days")),color=colors[i], label = labels[i])
    for k in range(len(profs[i].get("age"))):
        if profs[i].get("point_mass_index")[k] > 0:
            break
    axes.plot(profs[i].get("age")[k]/1e6,np.log10(profs[i].get("period_days")[k]),'o',color=colors[i], ms=6)

axes.text(4.6,0.75,"$M_{\\rm BH}=25,\\;a=0.35$", color = hexcols[2])
axes.text(4.6,0.45,"$M_{\\rm BH}=39,\\;a=0.77$", color = hexcols[8])
axes.text(4.6,0.26,"$M_{\\rm BH}=55,\\;a=1$", color = hexcols[5])

plt.text(3.5, 1.1, 'Widening due\nto mass loss',
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center', fontsize=8)
plt.text(15.5, 0.3, 'Widening due\nto mass transfer',
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center', fontsize=8)

axes.set_xlabel("age ${\\rm [Myr]}$")
axes.set_ylabel("$\\log P_{\\rm orb}\;{\\rm [days]}$")
axes.legend(loc="upper left", title = "$M_1=70M_\\odot,M_2=14M_\\odot$")

plt.savefig("../../images/period.pdf")

#        axarr[i].plot(profs_B[i].get("log_Teff")[3791:4273],profs_B[i].get("log_L")[3791:4273],color=hexcols[8],\
#                path_effects=[pe.Stroke(linewidth=7, foreground='k'), pe.Normal()], solid_capstyle='round',lw=6, zorder=-100)
