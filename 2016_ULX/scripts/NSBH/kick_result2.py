#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import *
import numpy as np
import sys
sys.path.insert(0, '../')
import kicks
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
          'lines.markersize': 4,
          'lines.linewidth': 3,
          'lines.antialiased': False,
          'path.simplify': False,
          'legend.handlelength':3,
          'figure.subplot.bottom':0.15,
          'figure.subplot.top':0.95,
          'figure.subplot.left':0.15,
          'figure.subplot.right':0.92}

hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',\
        '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466','#4477AA']

mpl.rcParams.update(params)

import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.set_cmap(cmaps.viridis)

data = np.loadtxt("kick.dataBH2")
orbit = kicks.post_kick_parameters_P(1.391271e+01, 3.635118e+01, 1.331844e+02, 3.635118e+01*0.9,0,0,0)

norm = mpl.colors.Normalize(vmin=0,vmax=16)

fig, axes= plt.subplots(1)
scatter = axes.scatter(data[:,3][np.logical_and(data[:,5]>0,True)],\
        data[:,4][np.logical_and(data[:,5]>0,True)], c=data[:,5][np.logical_and(data[:,5]>0,True)],\
        marker="o", s=10, linewidth='0', cmap = "viridis", norm=norm, rasterized = True)
axes.scatter(orbit[0], orbit[1], marker="s", s=40, linewidth='0')

for kickv in [50,100]:
    orbit = [kicks.post_kick_parameters_P(1.391271e+01, 3.635118e+01, 1.331844e+02, 3.635118e+01*0.9,kickv,theta,0) \
            for theta in np.linspace(0,math.pi,50)]
    orbit = np.array(orbit)
    axes.plot(orbit[:,0],orbit[:,1],"k-",linewidth=0.5)

#for kicktheta in [math.pi,math.pi-math.pi/4,math.pi-math.pi/3]:
#    orbit = [kicks.post_kick_parameters_P(12.7,48.8,136,43.9,kickv,kicktheta,0) \
#            for kickv in np.linspace(0,100,500)]
#    orbit = np.array(orbit)
#    axes.plot(orbit[:,0],orbit[:,1],color="0.5",linestyle="--",linewidth=0.5)

axes.text(15,0.17,"$v=50\\;{\\rm km\\;s^{-1}}$", rotation =25, fontsize = 7)
axes.text(17,0.36,"$v=100\\;{\\rm km\\;s^{-1}}$", rotation =20, fontsize = 7)

axes.text(1,0.92,"$P=14\\;{\\rm d}$", fontsize = 8)
axes.text(1,0.85,"$M_1=133M_\\odot$", fontsize = 8)
axes.text(1,0.78,"$M_2=36M_\\odot\\rightarrow 33M_\\odot$", fontsize = 8)

axes.text(24,0.92,"disrupt=$0\\%$", ha = "right", fontsize = 8)
axes.text(24,0.85,"merge=$0.014\\%$", ha = "right", fontsize = 8)

cbar = plt.colorbar(scatter)
cbar.set_label("Merger time ${\\rm[Gyr]}$")
axes.set_xlabel("$P_{\\rm orb}\\;{\\rm[days]}$")
axes.set_ylabel("eccentricity")
axes.set_xlim([0,25])
axes.set_ylim([0,1])
axes.legend(loc="best", scatterpoints=1)
#CS = plt.contour(lg_mass, vrot, He, levels=[0.3,0.5,0.7])
#plt.clabel(CS, inline=1, fontsize=10)

plt.savefig("kick_result2.pdf")
plt.clf()
plt.close(plt.gcf())
