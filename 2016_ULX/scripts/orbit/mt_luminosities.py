#!/usr/bin/env python
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
import matplotlib.patches as patches
from scipy.interpolate import griddata
import math
import scipy
from matplotlib import ticker
import sys
import os
import mmap
import itertools
import matplotlib as mpl
import matplotlib.gridspec as grd
from matplotlib.colors import LogNorm

params = {'backend': 'pdf',
          'figure.figsize': [4.3, 11],
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
          'figure.subplot.bottom':0.06,
          'figure.subplot.top':0.98,
          'figure.subplot.left':0.15,
          'figure.subplot.right':0.9}


mpl.rcParams.update(params)

clight = 3e10
cgrav = 6.67e-8
msun = 2e33
Lsun = 3.9e33

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

def clamp(val, minimum=0, maximum=255):
    if val < minimum:
        return minimum
    if val > maximum:
        return maximum
    return val

def colorscale(hexstr, scalefactor):
    """
    Scales a hex string by ``scalefactor``. Returns scaled hex string.

    To darken the color, use a float value between 0 and 1.
    To brighten the color, use a float value greater than 1.

    >>> colorscale("#DF3C3C", .5)
    #6F1E1E
    >>> colorscale("#52D24F", 1.6)
    #83FF7E
    >>> colorscale("#4F75D2", 1)
    #4F75D2
    """

    hexstr = hexstr.strip('#')

    if scalefactor < 0 or len(hexstr) != 6:
        return hexstr

    r, g, b = int(hexstr[:2], 16), int(hexstr[2:4], 16), int(hexstr[4:], 16)

    r = clamp(r * scalefactor)
    g = clamp(g * scalefactor)
    b = clamp(b * scalefactor)

    return "#%02x%02x%02x" % (r, g, b)

range_integral = 0.0194213341446518*(2.56+0.301)
limited_range = 0.004105169206845768*(0.5+0.301)

folders = ["-2.50","-3.00","-3.50","-4.00","-4.50","-5.00","-5.50","-6.00"]
folders = ["-2.500","-3.000","-3.500","-4.000","-4.500","-5.000","-5.500","-6.000"]
#folders = ["-4.500","-5.000"]
qratios = ["0.050","0.100","0.150","0.200","0.250","0.300","0.350","0.400","0.450","0.500","0.550","0.600"]
for folder in folders:
    luminosities = []
    spins = []
    superedd = []
    bh_masses = []
    donor_masses = []
    periods = []
    weights = []
    for qratio in qratios:
        file_name = "mt_data_"+folder+"_"+qratio+".dat"
        print(file_name)
        try:
            data = np.loadtxt(file_name, skiprows=2, unpack = True)
            if data.size == 0: continue
        except:
            print("cant open file, or no data", file_name)
            continue
        luminosities.append(data[11])
        spins.append(data[12])
        superedd.append(data[10]-data[9])
        bh_masses.append(data[6])
        donor_masses.append(data[7])
        periods.append(data[15])
        weights.append(data[4]*data[8])

    luminosities = np.log10(np.power(10,np.concatenate(luminosities))*Lsun)-39
    spins = np.concatenate(spins)
    superedd = np.concatenate(superedd)
    bh_masses = np.log10(np.concatenate(bh_masses))
    donor_masses = np.log10(np.concatenate(donor_masses))
    periods = np.log10(np.concatenate(periods))
    weights = np.concatenate(weights)

    sum_weights = np.sum(weights)
    weights = weights/sum_weights

    print("sources per SFR:", sum_weights*0.01/3)

    import colormaps as cmaps
    plt.register_cmap(name='viridis', cmap=cmaps.viridis)
    plt.set_cmap(cmaps.viridis)

    gs = grd.GridSpec(6, 3, height_ratios=[1,3,3,3,3,3], width_ratios=[20,5,1], wspace=0.1, hspace=0.1)

    axarr = [plt.subplot(gs[0]),plt.subplot(gs[4]),plt.subplot(gs[3]),plt.subplot(gs[5]),\
            plt.subplot(gs[7]), plt.subplot(gs[6]), plt.subplot(gs[8]),\
            plt.subplot(gs[10]), plt.subplot(gs[9]), plt.subplot(gs[11]),\
            plt.subplot(gs[13]), plt.subplot(gs[12]), plt.subplot(gs[14]),\
            plt.subplot(gs[16]), plt.subplot(gs[15]), plt.subplot(gs[17])]

    #histograms for mass-luminosity
    axarr[0].xaxis.set_major_locator(plt.NullLocator())
    axarr[0].yaxis.set_major_locator(plt.NullLocator())
    n, bins, patches = axarr[0].hist( bh_masses, weights=weights, color =hexcols[1],\
            rwidth = 1, bins=np.linspace(1.3-0.025,2.4+0.025,24),lw=0)
    axarr[0].set_xlim([1.3,2.4])
    axarr[0].set_xticklabels([])

    axarr[1].xaxis.set_major_locator(plt.NullLocator())
    axarr[1].yaxis.set_major_locator(plt.NullLocator())
    n, bins, patches = axarr[1].hist( luminosities, weights=weights, color =hexcols[1],\
            rwidth = 1, bins=np.linspace(39.2-0.33333-39,40.8+0.33333-39,26),lw=0, orientation="horizontal")
    axarr[1].set_ylim([39.2-39,40.7-39])
    axarr[1].set_yticklabels([])

    #density plot for mass-luminosity
    axarr[2].set_ylabel("$ \\log_{10} L_{\\rm acc}\\;\\mathrm{[10^{39}\;erg\;s^{-1}]}$", color="k")
    H, xedges, yedges = np.histogram2d(luminosities, bh_masses, weights = weights,\
            bins=(np.linspace(39.2-0.33333-39,40.8+0.33333-39,26), np.linspace(1.3-0.025,2.4+0.025,24)))
    X, Y = np.meshgrid(yedges, xedges)
    axarr[2].pcolormesh(X, Y, np.log10(H/np.max(H)), vmin=-4., vmax=0., linewidth=0,rasterized=True)
    axarr[2].set_xlim([1.3,2.4])
    axarr[2].set_ylim([39.2-39,40.7-39])
    axarr[2].set_xticklabels([])
    axarr[2].text(1.82,40.18-39,"$\\dot{L}_{\\rm Edd}(X=0.7)$", rotation=30, color="w",fontsize=10)
    axarr[2].text(1.73,40.53-39,"$\\dot{L}_{\\rm Edd}(X=0)$", rotation=30, color="w",fontsize=10)
    axarr[2].text(2.37,39.25-39,"$\\log_{10} Z="+folder[0:4]+"$", color="w",fontsize=15, ha = "right", va = "bottom")

    x = np.linspace(1,3,100)
    axarr[2].plot(x, np.log10(2.5e38/(1.7)) + x -39,"w--",linewidth=1)
    axarr[2].plot(x, np.log10(2.5e38/(1.0)) + x -39,"w--",linewidth=1)

    norm = norm=LogNorm(vmin=1e-4,vmax=1)
    mpl.colorbar.ColorbarBase(axarr[3], norm=norm, orientation='vertical')

    #histogram for spins
    axarr[4].xaxis.set_major_locator(plt.NullLocator())
    axarr[4].yaxis.set_major_locator(plt.NullLocator())
    n, bins, patches = axarr[4].hist( spins, weights=weights, color =hexcols[1],\
            rwidth = 1, bins=np.linspace(0,1,16),lw=0, orientation="horizontal")
    axarr[4].set_ylim([0,1])

    #density plot for spins
    axarr[5].set_ylabel("Black hole spin")
    H, xedges, yedges = np.histogram2d(spins, bh_masses, weights = weights,\
            bins=(np.linspace(0,1,16), np.linspace(1.3-0.025,2.4+0.025,24)))
    X, Y = np.meshgrid(yedges, xedges)
    print(np.max(H))
    axarr[5].pcolormesh(X, Y, np.log10(H/np.max(H)), vmin=-4., vmax=0., linewidth=0,rasterized=True)
    axarr[5].set_xlim([1.3,2.4])
    axarr[5].set_ylim([0,1])
    axarr[5].set_xticklabels([])
    mpl.colorbar.ColorbarBase(axarr[6], norm=norm, orientation='vertical')

    #histogram for superedd
    axarr[7].xaxis.set_major_locator(plt.NullLocator())
    axarr[7].yaxis.set_major_locator(plt.NullLocator())
    n, bins, patches = axarr[7].hist( superedd, weights=weights, color =hexcols[1],\
            rwidth = 1, bins=np.linspace(-0.25,2,17),lw=0, orientation="horizontal")
    axarr[7].set_ylim([-0.25,2])

    #density plot for superedd
    axarr[8].set_ylabel("$\\log_{10} \\dot{M}_{\\rm mt}/\\dot{M}_{\\rm Edd}$")
    H, xedges, yedges = np.histogram2d(superedd, bh_masses, weights = weights,\
            bins=(np.linspace(-0.25,2,17), np.linspace(1.3-0.025,2.4+0.025,24)))
    X, Y = np.meshgrid(yedges, xedges)
    print(np.max(H))
    axarr[8].pcolormesh(X, Y, np.log10(H/np.max(H)), vmin=-4., vmax=0., linewidth=0,rasterized=True)
    axarr[8].set_xlim([1.3,2.4])
    axarr[8].set_ylim([-0.25,2])
    axarr[8].set_xticklabels([])
    mpl.colorbar.ColorbarBase(axarr[9], norm=norm, orientation='vertical')

    #histogram for donor mass
    axarr[10].xaxis.set_major_locator(plt.NullLocator())
    axarr[10].yaxis.set_major_locator(plt.NullLocator())
    n, bins, patches = axarr[10].hist( donor_masses, weights=weights, color =hexcols[1],\
            rwidth = 1, bins=np.linspace(0,2.0,17),lw=0, orientation="horizontal")
    axarr[10].set_ylim([0,2.0])

    #density plot for donor mass
    axarr[11].set_ylabel("$\\log_{10}M_{\\rm donor}\;\mathrm{[M_\odot]}$")
    H, xedges, yedges = np.histogram2d(donor_masses, bh_masses, weights = weights,\
            bins=(np.linspace(0,2,17), np.linspace(1.3-0.025,2.4+0.025,24)))
    X, Y = np.meshgrid(yedges, xedges)
    print(np.max(H))
    axarr[11].pcolormesh(X, Y, np.log10(H/np.max(H)), vmin=-4., vmax=0., linewidth=0,rasterized=True)
    axarr[11].set_xlim([1.3,2.4])
    axarr[11].set_ylim([0,2])
    axarr[11].set_xticklabels([])
    mpl.colorbar.ColorbarBase(axarr[12], norm=norm, orientation='vertical')

    #histogram for donor mass
    axarr[13].xaxis.set_major_locator(plt.NullLocator())
    axarr[13].yaxis.set_major_locator(plt.NullLocator())
    n, bins, patches = axarr[13].hist( periods, weights=weights, color =hexcols[1],\
            rwidth = 1, bins=np.linspace(0,2.0,17),lw=0, orientation="horizontal")
    axarr[13].set_ylim([0,2.0])

    #density plot for donor mass
    axarr[14].set_ylabel("$\\log_{10}P_{\\rm orb}\;\mathrm{[days]}$")
    axarr[14].set_xlabel("$\\log_{10}M_{\mathrm{BH}}\;\mathrm{[M_\odot]}$")
    H, xedges, yedges = np.histogram2d(periods, bh_masses, weights = weights,\
            bins=(np.linspace(0,2,17), np.linspace(1.3-0.025,2.4+0.025,24)))
    X, Y = np.meshgrid(yedges, xedges)
    print(np.max(H))
    axarr[14].pcolormesh(X, Y, np.log10(H/np.max(H)), vmin=-4., vmax=0., linewidth=0,rasterized=True)
    axarr[14].set_xlim([1.3,2.4])
    axarr[14].set_ylim([0,2])
    mpl.colorbar.ColorbarBase(axarr[15], norm=norm, orientation='vertical')

    plt.savefig("plots/mt_luminosities_"+folder+".pdf", dpi=300, facecolor='w', edgecolor='w',
                    orientation='portrait', papertype=None, format=None,
                            transparent=False)
    plt.clf()
    plt.close(plt.gcf())
