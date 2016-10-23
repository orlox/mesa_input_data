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
          'figure.figsize': [4.3, 10],
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
          'figure.subplot.bottom':0.05,
          'figure.subplot.top':0.92,
          'figure.subplot.left':0.15,
          'figure.subplot.right':0.95}


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

lums_X = [\
        2.3,
        2.1,
        6.7,
        10.4,
        24.8,
        31.1,
        2.1,
        6.9,
        11.0,
        36.0,
        2.4,
        6.0,
        3.0,
        3.7,
        5.4,
        1.5,
        25.7,
        1.9,
        32.9,
        2.4,
        5.2,
        2.8,
        3.6,
        5.2,
        7.6,
        2.7,
        5.4,
        11.9,
        2.6,
        14.4,
        1.7,
        2.1,
        1.6,
        5.3,
        5.5,
        1.5,
        19.8,
        6.3,
        6.5,
        2.9,
        3.5,
        8.3,
        10.4,
        17.9,
        3.9,
        16.8,
        1.8,
        3.6,
        21.8,
        1.2,
        1.7,
        1.2,
        1.2,
        3.1,
        30.4,
        1.4,
        12.2,
        3.3,
        3.5,
        6.9,
        2.2,
        2.1,
        1.3,
        3.9,
        1.0,
        3.1,
        14.4,
        2.4,
        14.2,
        2.5,
        34.0,
        3.1,
        2.0,
        3.3,
        5.7,
        4.5,
        3.2,
        3.7,
        2.4]

lums_cts = [
        0.8,
        1.9,
        3.9,
        4.7,
        4.2,
        19.9,
        2.4,
        1.7,
        3.4,
        2.3,
        1.3,
        5.7,
        1.1,
        1.1,
        6.6,
        2.5,
        1.0,
        1.7,
        5.0,
        1.1,
        4.2,
        6.0,
        1.1,
        2.8,
        2.5,
        0.6,
        1.7,
        0.3,
        0.6,
        0.4,
        0.2,
        2.3,
        2.6,
        1.2,
        1.6,
        1.4,
        1.9,
        4.6,
        1.2,
        0.9,
        1.0,
        14.1,
        1.7,
        11.4,
        3.5,
        1.4,
        1.6,
        1.1,
        1.8,
        5.6,
        10.8,
        1.4,
        1.2,
        2.5,
        2.4,
        0.8,
        0.9,
        7.3,
        2.2,
        1.7,
        1.8,
        3.2,
        21.4,
        2.6,
        14.4,
        2.5,
        1.4,
        1.1,
        3.2,
        1.2,
        1.3,
        1.7,
        1.4,
        1.5,
        2.7,
        2.1,
        8.3,
        0.5,
        2.0,
        3.4,
        0.9,
        0.7,
        1.7,
        1.2,
        0.8,
        2.6,
        1.1,
        2.6,
        4.7,
        2.4,
        8.3,
        1.6,
        0.7,
        0.6,
        0.3,
        6.8,
        0.4,
        1.1,
        6.5,
        1.5,
        1.0,
        1.5,
        1.7]

Lvals_X = np.sort(lums_X)
nums_X = np.arange(len(Lvals_X),0,-1) / (30.3*1.7)

Lvals_cts = np.sort(lums_cts)
nums_cts = np.arange(len(Lvals_cts),0,-1) / (30.3*1.7)

Lgrimm = np.linspace(0.1,20)
Ygrimm = 5.4*(np.power(10*Lgrimm,-0.61)-210**(-0.61))

#fig, axes= plt.subplots(1)
gs = grd.GridSpec(8, 1, wspace=0.0, hspace=0.0)

folders = ["-2.500","-3.000","-3.500","-4.000","-4.500","-5.000","-5.500","-6.000"]
qratios = ["0.050","0.100","0.150","0.200","0.250","0.300","0.350","0.400","0.450","0.500","0.550","0.600"]

numsfr = [
        "0.6 ULX per $M_\odot~\\rm yr^{-1}$ of SFR",
        "1.9 ULX per $M_\odot~\\rm yr^{-1}$ of SFR",
        "2.2 ULX per $M_\odot~\\rm yr^{-1}$ of SFR",
        "2.3 ULX per $M_\odot~\\rm yr^{-1}$ of SFR",
        "1.5 ULX per $M_\odot~\\rm yr^{-1}$ of SFR",
        "1.1 ULX per $M_\odot~\\rm yr^{-1}$ of SFR",
        "0.56 ULX per $M_\odot~\\rm yr^{-1}$ of SFR",
        "0.11 ULX per $M_\odot~\\rm yr^{-1}$ of SFR"
        ]
numccsn = [
        "2.6 ULX per $10^4$ SNe",
        "1.9 ULX per $10^4$ SNe",
        "1.0 ULX per $10^4$ SNe",
        "1.0 ULX per $10^4$ SNe",
        "0.78 ULX per $10^4$ SNe",
        "0.66 ULX per $10^4$ SNe",
        "0.34 ULX per $10^4$ SNe",
        "0.067 ULX per $10^4$ SNe"
        ]
for i, folder in enumerate(folders):
    axes = plt.subplot(gs[i])
    luminosities = []
    spins = []
    superedd = []
    bh_masses = []
    donor_masses = []
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
        weights.append(data[4]*data[8])

    luminosities = np.power(10,np.concatenate(luminosities))*Lsun/1e39
    spins = np.concatenate(spins)
    superedd = np.concatenate(superedd)
    bh_masses = np.log10(np.concatenate(bh_masses))
    donor_masses = np.log10(np.concatenate(donor_masses))
    weights = np.concatenate(weights)

    sum_weights = np.sum(weights)
    #weights = weights/sum_weights

    print("sources per SFR:", sum_weights*0.01/3)

    Lbins = 10**np.linspace(-1,3,160)
    hist, bin_edges = np.histogram(luminosities, bins=Lbins, weights = weights*0.01/3)
    for k in range(len(hist)-2,-1,-1):
        hist[k] = hist[k] + hist[k+1]
    left,right = bin_edges[:-1],bin_edges[1:]
    X_Z = np.array([left,right]).T.flatten()
    Y_Z = np.array([hist,hist]).T.flatten()
    scale_array = np.power(10,\
            np.minimum(np.maximum(superedd,0),np.log10(3)))
    hist2, bin_edges = np.histogram(luminosities*scale_array, bins=Lbins, weights = weights*0.01/3)
    for k in range(len(hist)-2,-1,-1):
        hist2[k] = hist2[k] + hist2[k+1]
    left,right = bin_edges[:-1],bin_edges[1:]
    X_Z2 = np.array([left,right]).T.flatten()
    Y_Z2 = np.array([hist2,hist2]).T.flatten()

    #axes.plot(X_Z,Y_Z, color=hexcols[6], label="This work")
    #axes.plot(X_Z2,Y_Z2, color="0.8")
    axes.fill_between(X_Z,1e-5*np.ones(len(Y_Z)),Y_Z, color=hexcols[6], alpha = 0.4, label="This work", linewidth = 3)
    axes.fill_between(X_Z2,1e-5*np.ones(len(Y_Z2)),Y_Z2, linestyle=":", color=hexcols[6], alpha = 0.2, label="This work ($3\\times \dot{M}_{\\rm Edd}$)", linewidth = 3)
    axes.plot(Lvals_X, nums_X, color=hexcols[1], ls="steps", label="Swartz et al. (2011), $L_{\\rm X}$")
    axes.plot(Lvals_cts, nums_cts, color =hexcols[3], ls="steps", label="Swartz et al. (2011), $L_{\\rm cnt}$")
    axes.plot(Lgrimm, Ygrimm, color =hexcols[0], label="Grimm et al. (2003)", linestyle="--")

    if i==0:
        axes.legend(loc='upper left', bbox_to_anchor=(-0.15, 1.7),
                ncol=2, scatterpoints=1,prop={'size':9})

    axes.set_yscale("log")
    axes.set_xscale("log")
    axes.set_ylim(0.002,4)
    axes.set_xlim(0.5,200)
    axes.set_ylabel("$N(>L)$")
    axes.text(40,0.7,"$\\log Z="+folder[:4]+"$")
    axes.text(0.6,0.012,numsfr[i], fontsize = 9)
    axes.text(0.6,0.004,numccsn[i], fontsize = 9)
    if i!= 7:
        axes.set_xticklabels([])
axes.set_xlabel("$L\\;[10^{39}\\;\\rm erg\\; s^{-1}]$")

plt.savefig("plots/lum_dist_cum.pdf", dpi=300, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                        transparent=False)
