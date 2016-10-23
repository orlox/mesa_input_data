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

data = np.loadtxt("kick.dataNS")
orbit = kicks.post_kick_parameters_P(2.348553e+01, 5.486435e+00, 5.504611e+01, 1.4,0,0,0)

norm = mpl.colors.Normalize(vmin=0,vmax=16)

fig, axes= plt.subplots(1)
scatter = axes.scatter(data[:,3][np.logical_and(data[:,5]>0,True)],\
        data[:,4][np.logical_and(data[:,5]>0,True)], c=data[:,5][np.logical_and(data[:,5]>0,True)],\
        marker="o", s=5, linewidth='0', cmap = "viridis", norm=norm, rasterized = True)
axes.scatter(orbit[0], orbit[1], marker="s", s=20, linewidth='0')

print("fraction below 1Gyr", \
        np.sum(data[:,6][np.logical_and(data[:,5]>0,data[:,5]<1)])\
        /np.sum(data[:,6][np.logical_and(data[:,5]>0,data[:,5]<13.8)]),
        np.sum(data[:,6][np.logical_and(data[:,5]>0,True)]),
        np.sum(data[:,6][np.logical_and(data[:,5]>0,data[:,5]<13.8)]))

print("fraction below 2Gyr", \
        np.sum(data[:,6][np.logical_and(data[:,5]>0,data[:,5]<2)])\
        /np.sum(data[:,6][np.logical_and(data[:,5]>0,data[:,5]<13.8)]),
        np.sum(data[:,6][np.logical_and(data[:,5]>0,True)]),
        np.sum(data[:,6][np.logical_and(data[:,5]>0,data[:,5]<13.8)]))

for kickv in [100,200,300]:#,400]:
    orbit = [kicks.post_kick_parameters_P(2.348553e+01, 5.486435e+00, 5.504611e+01, 1.4,kickv,theta,0) \
            for theta in np.linspace(0,math.pi,50)]
    orbit = np.array(orbit)
    axes.plot(orbit[:,0],orbit[:,1],"k-",linewidth=0.5)

for kicktheta in [math.pi,math.pi-math.pi/4,math.pi-math.pi/3]:
    orbit = [kicks.post_kick_parameters_P(2.348553e+01, 5.486435e+00, 5.504611e+01, 1.4,kickv,kicktheta,0) \
            for kickv in np.linspace(0,1000,500)]
    orbit = np.array(orbit)
    axes.plot(orbit[:,0],orbit[:,1],color="0.5",linestyle="--",linewidth=0.5)

axes.text(30,0.38,"$v=100\\;{\\rm km\\;s^{-1}}$", rotation =15, fontsize = 7)
axes.text(30,0.61,"$v=200\\;{\\rm km\\;s^{-1}}$", rotation =3, fontsize = 7)
axes.text(30,0.83,"$v=300\\;{\\rm km\\;s^{-1}}$", rotation =-3, fontsize = 7)

axes.text(8,0.75,"$\\theta=\\pi$", rotation =100, fontsize = 7)
axes.text(11.5,0.7,"$\\theta=3\\pi/4$", rotation =100, fontsize = 7)
axes.text(16.5,0.65,"$\\theta=2\\pi/3$", rotation =90, fontsize = 7)

axes.text(1,0.2,"$P_{\\rm f}=23\\;{\\rm d}$", fontsize = 8)
axes.text(1,0.13,"$M_{1,\\rm f}=55M_\\odot$", fontsize = 8)
axes.text(1,0.06,"$M_{2,\\rm f}=5.5M_\\odot\\rightarrow 1.4M_\\odot$", fontsize = 8)

axes.text(49,0.13,"disrupt=$67\\%$", ha = "right", fontsize = 8)
axes.text(49,0.06,"merge=$3.4\\%$", ha = "right", fontsize = 8)

cbar = plt.colorbar(scatter)
cbar.set_label("Merger time ${\\rm[Gyr]}$")
axes.set_xlabel("$P_{\\rm orb}\\;{\\rm[days]}$")
axes.set_ylabel("eccentricity")
axes.set_xlim([0,50])
axes.set_ylim([0,1])
axes.legend(loc="best", scatterpoints=1)
#CS = plt.contour(lg_mass, vrot, He, levels=[0.3,0.5,0.7])
#plt.clabel(CS, inline=1, fontsize=10)

plt.savefig("kick_result.pdf")
plt.clf()
plt.close(plt.gcf())
