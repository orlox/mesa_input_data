#!/usr/bin/env python
from mesa_data import *
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import matplotlib.patheffects as pe
import matplotlib.gridspec as grd
params = {'backend': 'pdf',
          'figure.figsize': [4.3, 6],
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
          'figure.subplot.top':0.85,
          'figure.subplot.left':0.15,
          'figure.subplot.right':0.92}

hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',\
        '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466','#4477AA']

mpl.rcParams.update(params)

Lsun = 3.9e33

plt.figure()
gs = grd.GridSpec(3,3,height_ratios=[1,1,1], wspace=0.05, hspace=0.05)
axes0 = plt.subplot(gs[6])
axes0b = plt.subplot(gs[7])
axes0c = plt.subplot(gs[8])
axes1 = plt.subplot(gs[3])
axes1b = plt.subplot(gs[4])
axes1c = plt.subplot(gs[5])
axes2 = plt.subplot(gs[0])
axes2b = plt.subplot(gs[1])
axes2c = plt.subplot(gs[2])

profs = [Mesa_Data("results_normal/binary_history.data"), Mesa_Data("results_3times/binary_history.data"), Mesa_Data("results_10times/binary_history.data")]
#colors = [hexcols[2], hexcols[8], hexcols[5]]
#labels = ["$\\log Z=-2.5, P_{\\rm orb,i}=1.4{\\rm d}$",\
#        "$\\log Z=-3.0, P_{\\rm orb,i}=1.2{\\rm d}$",\
#        "$\\log Z=-3.5, P_{\\rm orb,i}=1.1{\\rm d}$"]
axes2.plot(profs[2].get("age")/1e6,np.power(10,profs[2].get("lg_mstar_dot_1"))/np.power(10,profs[2].get("lg_mstar_dot_2")),color=hexcols[2], label = "$10\\times\\dot{M}_{\\rm Edd}$ limited")
axes2.plot(profs[1].get("age")/1e6,np.power(10,profs[1].get("lg_mstar_dot_1"))/np.power(10,profs[1].get("lg_mstar_dot_2")),color=hexcols[5], label = "$3\\times\\dot{M}_{\\rm Edd}$ limited")
axes2.plot(profs[0].get("age")/1e6,np.power(10,profs[0].get("lg_mstar_dot_1"))/np.power(10,profs[0].get("lg_mstar_dot_2")),color=hexcols[7], label = "$\\dot{M}_{\\rm Edd}$ limited")
axes2b.plot(profs[2].get("age")/1e6,np.power(10,profs[2].get("lg_mstar_dot_1"))/np.power(10,profs[2].get("lg_mstar_dot_2")),color=hexcols[2], label = "$10\\times\\dot{M}_{\\rm Edd}$ limited")
axes2b.plot(profs[1].get("age")/1e6,np.power(10,profs[1].get("lg_mstar_dot_1"))/np.power(10,profs[1].get("lg_mstar_dot_2")),color=hexcols[5], label = "$10\\times\\dot{M}_{\\rm Edd}$ limited")
axes2b.plot(profs[0].get("age")/1e6,np.power(10,profs[0].get("lg_mstar_dot_1"))/np.power(10,profs[0].get("lg_mstar_dot_2")),color=hexcols[7], label = "$10\\times\\dot{M}_{\\rm Edd}$ limited")
axes2c.plot(profs[2].get("age")/1e6,np.power(10,profs[2].get("lg_mstar_dot_1"))/np.power(10,profs[2].get("lg_mstar_dot_2")),color=hexcols[2], label = "$10\\times\\dot{M}_{\\rm Edd}$ limited")
axes2c.plot(profs[1].get("age")/1e6,np.power(10,profs[1].get("lg_mstar_dot_1"))/np.power(10,profs[1].get("lg_mstar_dot_2")),color=hexcols[5], label = "$10\\times\\dot{M}_{\\rm Edd}$ limited")
axes2c.plot(profs[0].get("age")/1e6,np.power(10,profs[0].get("lg_mstar_dot_1"))/np.power(10,profs[0].get("lg_mstar_dot_2")),color=hexcols[7], label = "$10\\times\\dot{M}_{\\rm Edd}$ limited")

axes1.plot(profs[2].get("age")/1e6,np.log10(profs[2].get("period_days")),color=hexcols[2], label = "$10\\times\\dot{M}_{\\rm Edd}$ limited")
axes1.plot(profs[1].get("age")/1e6,np.log10(profs[1].get("period_days")),color=hexcols[5], label = "$3\\times\\dot{M}_{\\rm Edd}$ limited")
axes1.plot(profs[0].get("age")/1e6,np.log10(profs[0].get("period_days")),color=hexcols[7], label = "$\\dot{M}_{\\rm Edd}$ limited")
axes1b.plot(profs[2].get("age")/1e6,np.log10(profs[2].get("period_days")),color=hexcols[2], label = "$10\\times\\dot{M}_{\\rm Edd}$ limited")
axes1b.plot(profs[1].get("age")/1e6,np.log10(profs[1].get("period_days")),color=hexcols[5], label = "$3\\times\\dot{M}_{\\rm Edd}$ limited")
axes1b.plot(profs[0].get("age")/1e6,np.log10(profs[0].get("period_days")),color=hexcols[7], label = "$\\dot{M}_{\\rm Edd}$ limited")
axes1c.plot(profs[2].get("age")/1e6,np.log10(profs[2].get("period_days")),color=hexcols[2], label = "$10\\times\\dot{M}_{\\rm Edd}$ limited")
axes1c.plot(profs[1].get("age")/1e6,np.log10(profs[1].get("period_days")),color=hexcols[5], label = "$3\\times\\dot{M}_{\\rm Edd}$ limited")
axes1c.plot(profs[0].get("age")/1e6,np.log10(profs[0].get("period_days")),color=hexcols[7], label = "$\\dot{M}_{\\rm Edd}$ limited")

axes0.plot(profs[2].get("age")/1e6,profs[2].get("lg_accretion_luminosity")-39+np.log10(Lsun),\
        color=hexcols[2], label = "$10 \\times \\dot{M}_{\\rm Edd}$")
axes0.plot(profs[0].get("age")/1e6,
        np.log10(Lsun*np.power(10,profs[0].get("lg_accretion_luminosity"))\
            *np.minimum(np.maximum(
                np.ones(\
                        len(profs[0].get("age"))),\
                        np.power(10,profs[0].get("lg_mtransfer_rate")-profs[0].get("lg_mdot_edd"))\
                ), 10)\
        )-39,\
        color="0.3", ls=":", label = "$10 \\times \\dot{M}_{\\rm Edd}$, from $\\dot{M}_{\\rm Edd}$", lw=1.5)

axes0.plot(profs[1].get("age")/1e6,profs[1].get("lg_accretion_luminosity")-39+np.log10(Lsun),\
        color=hexcols[5], label = "$3\\times \\dot{M}_{\\rm Edd}$")
axes0.plot(profs[0].get("age")/1e6,
        np.log10(Lsun*np.power(10,profs[0].get("lg_accretion_luminosity"))\
            *np.minimum(np.maximum(
                np.ones(\
                        len(profs[0].get("age"))),\
                        np.power(10,profs[0].get("lg_mtransfer_rate")-profs[0].get("lg_mdot_edd"))\
                ), 3)\
        )-39,\
        color="k", ls="--", label = "$3\\times \\dot{M}_{\\rm Edd}$, from $\\dot{M}_{\\rm Edd}$", lw=1.5)

axes0.plot(profs[0].get("age")/1e6,profs[0].get("lg_accretion_luminosity")-39+np.log10(Lsun),\
        color=hexcols[7], label = "$10 \\times \\dot{M}_{\\rm Edd}$")

axes0b.plot(profs[2].get("age")/1e6,profs[2].get("lg_accretion_luminosity")-39+np.log10(Lsun),\
        color=hexcols[2], label = "$10 \\times \\dot{M}_{\\rm Edd}$")
axes0b.plot(profs[0].get("age")/1e6,
        np.log10(Lsun*np.power(10,profs[0].get("lg_accretion_luminosity"))\
            *np.minimum(np.maximum(
                np.ones(\
                        len(profs[0].get("age"))),\
                        np.power(10,profs[0].get("lg_mtransfer_rate")-profs[0].get("lg_mdot_edd"))\
                ), 10)\
        )-39,\
        color="0.3", ls=":", label = "$10 \\times \\dot{M}_{\\rm Edd}$, from $\\dot{M}_{\\rm Edd}$", lw=1.5)

axes0b.plot(profs[1].get("age")/1e6,profs[1].get("lg_accretion_luminosity")-39+np.log10(Lsun),\
        color=hexcols[5], label = "$3\\times \\dot{M}_{\\rm Edd}$")
axes0b.plot(profs[0].get("age")/1e6,
        np.log10(Lsun*np.power(10,profs[0].get("lg_accretion_luminosity"))\
            *np.minimum(np.maximum(
                np.ones(\
                        len(profs[0].get("age"))),\
                        np.power(10,profs[0].get("lg_mtransfer_rate")-profs[0].get("lg_mdot_edd"))\
                ), 3)\
        )-39,\
        color="k", ls="--", label = "$3\\times \\dot{M}_{\\rm Edd}$, from $\\dot{M}_{\\rm Edd}$", lw=1.5)
axes0b.plot(profs[0].get("age")/1e6,profs[0].get("lg_accretion_luminosity")-39+np.log10(Lsun),\
        color=hexcols[7], label = "$10 \\times \\dot{M}_{\\rm Edd}$")

axes0c.plot(profs[2].get("age")/1e6,profs[2].get("lg_accretion_luminosity")-39+np.log10(Lsun),\
        color=hexcols[2], label = "$10 \\times \\dot{M}_{\\rm Edd}$")
axes0c.plot(profs[0].get("age")/1e6,
        np.log10(Lsun*np.power(10,profs[0].get("lg_accretion_luminosity"))\
            *np.minimum(np.maximum(
                np.ones(\
                        len(profs[0].get("age"))),\
                        np.power(10,profs[0].get("lg_mtransfer_rate")-profs[0].get("lg_mdot_edd"))\
                ), 10)\
        )-39,\
        color="0.3", ls=":", label = "$10 \\times \\dot{M}_{\\rm Edd}$, from $\\dot{M}_{\\rm Edd}$", lw=1.5)

axes0c.plot(profs[1].get("age")/1e6,profs[1].get("lg_accretion_luminosity")-39+np.log10(Lsun),\
        color=hexcols[5], label = "$3\\times \\dot{M}_{\\rm Edd}$")
axes0c.plot(profs[0].get("age")/1e6,
        np.log10(Lsun*np.power(10,profs[0].get("lg_accretion_luminosity"))\
            *np.minimum(np.maximum(
                np.ones(\
                        len(profs[0].get("age"))),\
                        np.power(10,profs[0].get("lg_mtransfer_rate")-profs[0].get("lg_mdot_edd"))\
                ), 3)\
        )-39,\
        color="k", ls="--", label = "$3\\times \\dot{M}_{\\rm Edd}$, from $\\dot{M}_{\\rm Edd}$", lw=1.5)
axes0c.plot(profs[0].get("age")/1e6,profs[0].get("lg_accretion_luminosity")-39+np.log10(Lsun),\
        color=hexcols[7], label = "$10 \\times \\dot{M}_{\\rm Edd}$")

#axes1.text(13.7, -5, 'Case A',
#         horizontalalignment='center',
#         verticalalignment='top',
#         multialignment='center', fontsize=7)
#axes1.text(15.5, -3, 'Case AB/B',
#         horizontalalignment='center',
#         verticalalignment='top',
#         multialignment='center', fontsize=7)
#axes1.text(16.2, -3.9, 'Case ABB/BB',
#         horizontalalignment='center',
#         verticalalignment='top',
#         multialignment='center', fontsize=7)
#
#axes0.text(13.2, 39.8-39, '$\dot{M}_{\\rm edd}$ limited',
#         horizontalalignment='center',
#         verticalalignment='top',
#         multialignment='center', fontsize=7)
#axes0.text(13.3, 41-39, '$\dot{M}_{\\rm edd}$ limit ignored',
#         horizontalalignment='center',
#         verticalalignment='top',
#         multialignment='center', fontsize=7, rotation = 8)
xlim1 = [13.6,15.3]
xlim2 = [15.32,15.41]
xlim3 = [16.29,16.35]

axes2.set_xlim(xlim1)
axes2b.set_xlim(xlim2)
axes2c.set_xlim(xlim3)
axes2.set_ylim([-0.05,1.05])
axes2b.set_ylim([-0.05,1.05])
axes2c.set_ylim([-0.05,1.05])
axes1.set_xlim(xlim1)
axes1b.set_xlim(xlim2)
axes1c.set_xlim(xlim3)
axes1.set_ylim([0.3,1.7])
axes1b.set_ylim([0.3,1.7])
axes1c.set_ylim([0.3,1.7])
#axes1.set_ylim([-7,-2.2])
#axes1.yaxis.set_ticklabels([-7,-6,-5,-4,-3])
axes0.set_xlim(xlim1)
axes0b.set_xlim(xlim2)
axes0c.set_xlim(xlim3)
axes0.set_ylim([0,2.2])
axes0b.set_ylim([0,2.2])
axes0c.set_ylim([0,2.2])
axes1.xaxis.set_ticklabels([])
axes1b.xaxis.set_ticklabels([])
axes1c.xaxis.set_ticklabels([])
axes1b.yaxis.set_ticklabels([])
axes1c.yaxis.set_ticklabels([])
axes2.xaxis.set_ticklabels([])
axes2b.xaxis.set_ticklabels([])
axes2c.xaxis.set_ticklabels([])
axes2b.yaxis.set_ticklabels([])
axes2c.yaxis.set_ticklabels([])
axes0b.yaxis.set_ticklabels([])
axes0c.yaxis.set_ticklabels([])

axes0b.set_xlabel("age ${\\rm [Myr]}$")
axes2.set_ylabel("$-\dot{M}_{\\rm BH}/\dot{M}_{\\rm donor}$")
axes1.set_ylabel("$\\log P_{\\rm orb}\;{\\rm[days]}$")
axes0.set_ylabel("$\\log L_{\\rm acc}\;{[10^{39}\\rm erg\;s^{-1}]}$")

axes2.plot([],[],
        color="0.3", ls=":", label = "$10 \\times \\dot{M}_{\\rm Edd}$, from $\\dot{M}_{\\rm Edd}$")
axes2.plot([],[],
        color="k", ls="--", label = "$3 \\times \\dot{M}_{\\rm Edd}$, from $\\dot{M}_{\\rm Edd}$")
axes2.legend(loc='upper left', bbox_to_anchor=(-0.05, 1.6),
        ncol=2, scatterpoints=1,prop={'size':8},
        title = "$M_1=70M_\\odot,M_2=14M_\\odot,\\log Z=-3.0$")

plt.savefig("mdotedd_check2.pdf")

#        axarr[i].plot(profs_B[i].get("log_Teff")[3791:4273],profs_B[i].get("log_L")[3791:4273],color=hexcols[8],\
#                path_effects=[pe.Stroke(linewidth=7, foreground='k'), pe.Normal()], solid_capstyle='round',lw=6, zorder=-100)
