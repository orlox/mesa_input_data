#!/usr/bin/env python
from mesa_data import *
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import matplotlib.patheffects as pe
import matplotlib.gridspec as grd
params = {'backend': 'pdf',
          'figure.figsize': [4.3, 5.5],
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
          'figure.subplot.bottom':0.075,
          'figure.subplot.top':0.95,
          'figure.subplot.left':0.15,
          'figure.subplot.right':0.92}

hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',\
        '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466','#4477AA']

mpl.rcParams.update(params)

Lsun = 3.9e33

plt.figure()
gs = grd.GridSpec(2,1, wspace=0.05, hspace=0.05)
axes0 = plt.subplot(gs[1])
axes1 = plt.subplot(gs[0])

profs = [Mesa_Data("binary_history_250.data"), Mesa_Data("binary_history_300.data"), Mesa_Data("binary_history_350.data")]
colors = [hexcols[2], hexcols[8], hexcols[5]]
labels = ["$\\log Z=-2.5, P_{\\rm orb,i}=1.4{\\rm d}$",\
        "$\\log Z=-3.0, P_{\\rm orb,i}=1.2{\\rm d}$",\
        "$\\log Z=-3.5, P_{\\rm orb,i}=1.1{\\rm d}$"]
for i in [0,1,2]:
    axes1.plot(profs[i].get("age")/1e6,profs[i].get("lg_mtransfer_rate"),color=colors[i], label = labels[i], lw=1.5)
for i in [0,1,2]:
    axes1.plot(profs[i].get("age")/1e6,profs[i].get("lg_mdot_edd"),"--",lw=1,color=colors[i])

for i in [0,1,2]:
    axes0.fill_between(profs[i].get("age")/1e6,profs[i].get("lg_accretion_luminosity")+np.log10(Lsun)-39,\
            np.log10(Lsun*np.power(10,profs[i].get("lg_accretion_luminosity"))\
                *np.maximum(
                    np.ones(\
                            len(profs[i].get("age"))),\
                            np.power(10,profs[i].get("lg_mtransfer_rate")-profs[i].get("lg_mdot_edd"))\
                    )\
            )-39,\
            color=colors[i], label = labels[i], alpha=0.5)

axes1.plot([],[],ls="--",color="k", lw=1,label="$\\dot{M}_{\\rm Edd}$")

axes1.text(13.7, -5, 'Case A',
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center', fontsize=7)
axes1.text(15.5, -3, 'Case AB/B',
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center', fontsize=7)
axes1.text(16.2, -3.9, 'Case ABB/BB',
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center', fontsize=7)

axes0.text(13.2, 39.8-39, '$\dot{M}_{\\rm edd}$ limited',
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center', fontsize=7)
axes0.text(13.3, 41-39, '$\dot{M}_{\\rm edd}$ limit ignored',
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center', fontsize=7, rotation = 8)

axes1.set_xlim([12.5,16.7])
axes1.set_ylim([-7,-2.2])
axes1.yaxis.set_ticklabels([-7,-6,-5,-4,-3])
axes0.set_xlim([12.5,16.7])
axes0.set_ylim([0,4])
axes1.xaxis.set_ticklabels([])

axes0.set_xlabel("age ${\\rm [Myr]}$")
axes1.set_ylabel("$\\log\dot{M}_{\\rm mt}\;{[M_\\odot\;\\rm yr^{-1}]}$")
axes0.set_ylabel("$\\log L_{\\rm acc}\;{[10^{39}\\rm erg\;s^{-1}]}$")
axes0.legend(loc="upper left", title = "$M_1=70M_\\odot,M_2=14M_\\odot$")
axes1.legend(loc="upper left")

plt.savefig("../../images/masstransfer.pdf")

#        axarr[i].plot(profs_B[i].get("log_Teff")[3791:4273],profs_B[i].get("log_L")[3791:4273],color=hexcols[8],\
#                path_effects=[pe.Stroke(linewidth=7, foreground='k'), pe.Normal()], solid_capstyle='round',lw=6, zorder=-100)
