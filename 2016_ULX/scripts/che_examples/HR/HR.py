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

profs_A = [Mesa_Data("history005_A.data"), Mesa_Data("history020_A.data"), Mesa_Data("history060_A.data")]
colors = [hexcols[1], hexcols[5], hexcols[3]]
labels = ["$M_1=70M_\odot, q=0.05, P_{\\rm i}=0.8{\\rm d}$",\
        "$M_1=70M_\odot, q=0.2, P_{\\rm i}=1.1{\\rm d}$",\
        "$M_1=70M_\odot, q=0.6, P_{\\rm i}=1.2{\\rm d}$"]
for i in [1,2]:
    for k in range(len(profs_A[i].get("log_Teff"))):
        if profs_A[i].get("star_age")[k] > 50000:
            break
    axes.plot(profs_A[i].get("log_Teff")[k:],profs_A[i].get("log_L")[k:],color=colors[i], label = labels[i])
    axes.plot(profs_A[i].get("log_Teff")[-1],profs_A[i].get("log_L")[-1],'o',color=colors[i], ms=6)
    if i == 1:
        for j, centerh1 in enumerate(profs_A[i].get("center_h1")):
            if centerh1 < 1e-6:
                axes.plot(profs_A[i].get("log_Teff")[j],profs_A[i].get("log_L")[j],'o',color=colors[i], ms=6)
                break

axes.text(5.0,5.95,"Secondary RLOF", color = hexcols[3])
axes.text(5.29,6.45,"BH formation", color = hexcols[5])
axes.text(4.995,6.3,"TAMS", color = hexcols[5])
axes.text(4.78,5.76,"ZAMS")
axes.set_xlim([5.3,4.6])
axes.set_ylim([5.72,6.5])
axes.set_xlabel("$\\log~T_\\mathrm{eff}\\;\\rm[K]$")
axes.set_ylabel("$\\log~L\\;\\rm[L_\\odot]$")
axes.legend(loc="lower left")

text(0.95, 0.95,'Primary',
     horizontalalignment='right',
     verticalalignment='top',
     transform = axes.transAxes, fontsize=15)
plt.savefig("../../images/HR_primary.pdf")

plt.figure()
fig, axes= plt.subplots(1)

profs_B = [Mesa_Data("history005_B.data"), Mesa_Data("history020_B.data"), Mesa_Data("history060_B.data")]
labels = ["$M_1=70M_\odot, q=0.05, P_{\\rm i}=0.8{\\rm d}$",\
        "$M_1=70M_\odot, q=0.2, P_{\\rm i}=1.1{\\rm d}$",\
        "$M_1=70M_\odot, q=0.6, P_{\\rm i}=1.2{\\rm d}$"]
agelims = [500000,50000,50000]

for i in [1,2]:
    for k in range(len(profs_B[i].get("log_Teff"))):
        if profs_B[i].get("star_age")[k] > agelims[i]:
            break
    axes.plot(profs_B[i].get("log_Teff")[k:],profs_B[i].get("log_L")[k:],color=colors[i])
    axes.plot(profs_B[i].get("log_Teff")[-1],profs_B[i].get("log_L")[-1],'o',color=colors[i], ms=6)
    if i == 1:
        for k, centerhe4 in enumerate(profs_A[i].get("center_he4")):
            if centerhe4 < 1e-3:
                axes.plot(profs_B[i].get("log_Teff")[k],profs_B[i].get("log_L")[k],'o',color=colors[i], ms=6)
                break
        for k, centerh1 in enumerate(profs_B[i].get("center_h1")):
            if centerh1 < 1e-6:
                axes.plot(profs_B[i].get("log_Teff")[k],profs_B[i].get("log_L")[k],'o',color=colors[i], ms=6)
                break
        for k, centerhe4 in enumerate(profs_B[i].get("center_he4")):
            if centerhe4 < 1e-6:
                axes.plot(profs_B[i].get("log_Teff")[k],profs_B[i].get("log_L")[k],'o',color=colors[i], ms=6)
                break
        axes.plot(profs_B[i].get("log_Teff")[3791:4273],profs_B[i].get("log_L")[3791:4273],color=hexcols[8],\
                path_effects=[pe.Stroke(linewidth=7, foreground='k'), pe.Normal()], solid_capstyle='round',lw=6, zorder=-100)
        axes.plot(profs_B[i].get("log_Teff")[4471:4600],profs_B[i].get("log_L")[4471:4600],color=hexcols[8],\
                path_effects=[pe.Stroke(linewidth=7, foreground='k'), pe.Normal()], solid_capstyle='round',lw=6, zorder=-100)
        axes.plot(profs_B[i].get("log_Teff")[5891:6200],profs_B[i].get("log_L")[5891:6200],color=hexcols[8],\
                path_effects=[pe.Stroke(linewidth=7, foreground='k'), pe.Normal()], solid_capstyle='round',lw=6, zorder=-100)


axes.text(4.71,5.35,"ZAMS", color = hexcols[3],ha='right',va='center')
axes.text(4.67,5.5,"RLOF", color = hexcols[3],ha='left',va='center')
axes.text(4.65,4.2,"ZAMS", color = hexcols[5])
axes.text(4.53,4.3,"Primary forms BH", color = hexcols[5])
axes.text(4.44,4.44,"Case A", color = hexcols[8], rotation=-7)
axes.text(4.27,4.85,"Case AB", color = hexcols[8], rotation=80)
axes.text(4.28,5.0,"Case ABB", color = hexcols[8])
axes.text(4.4,4.65,"TAMS", color = hexcols[5])
axes.text(4.68,5.05,"He depletion", color = hexcols[5])
axes.text(4.37,5.1,"C depletion", color = hexcols[5])
#axes.set_ylim([4.1,5.2])
#axes.set_xlim([4.8,4.15])
axes.set_ylim([4.15,5.57])
axes.set_xlim([4.8,4.15])
axes.set_xlabel("$\\log~T_\\mathrm{eff}\;\\rm[K]$")
axes.set_ylabel("$\\log~L\\;\\rm[L_\\odot]$")
axes.legend(loc="lower left")

plt.sca(axes)
text(0.95, 0.95,'Secondary',
     horizontalalignment='right',
     verticalalignment='top',
     transform = axes.transAxes, fontsize=15)
plt.savefig("../../images/HR_secondary.pdf")
