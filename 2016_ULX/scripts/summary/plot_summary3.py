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
import matplotlib.gridspec as grd

params = {'backend': 'pdf',
          'figure.figsize': [4.3, 8],
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
          'figure.subplot.bottom':0.07,
          'figure.subplot.top':0.93,
          'figure.subplot.left':0.15,
          'figure.subplot.right':0.95}

mpl.rcParams.update(params)

clight = 3e10
cgrav = 6.67e-8
msun = 2e33

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

def plot_summary(ax,Pwidth, logMwidth, metallicity, qratio, file_name ='.', kick_data = None):
    """
    Plots a summary of a grid of binary models for fixed donor mass:
        mass_donor: string with the donor mass
        masses_accretor: array accretor masses, must be sorted in increasing order
        periods: array of strings with periods, must be sorted in increasing order
        hcols: dictionary with column numbers (starting from zero) of history file
        folder: base folder where models are located
    
    The folder with the data is assumed to be in:
    folder+mass_donor+"_"+mass_accretor+"_"+period
    
    Output is stored in:
    folder+"/summary_"+mass_donor+".pdf"
    
    """

    #data = np.genfromtxt(folder, dtype=None, names=['log10(M_1i)(Msun)','P_i(days)','result', 'had_contact'], skip_header=15)
    data = np.genfromtxt(file_name, dtype=None, names=True, skip_header=26)
    print(file_name)

    for datum in data:
        #try:
        #   if float(datum[11])>0:
        #       hatch = "xxxx"
        #except ValueError:
        #    hatch = ""
        hatch = ""

        if datum[4].decode('utf-8')=="ZAMS_L2_overflow":
            state = "black"
        elif datum[4].decode('utf-8')=="ZAMS_RLOF":
            state = "0.5"
        elif datum[4].decode('utf-8') == "off_CHE":
            state = hexcols[6]
        elif datum[4].decode('utf-8') == "MT_P->S":
            state = hexcols[2]
        elif datum[4].decode('utf-8') == "PISN":
            state = hexcols[12]
        elif datum[4][:-2].decode('utf-8') == "MT_to_pm_caseB":
            state = hexcols[8]
        elif datum[4][:-2].decode('utf-8') == "MT_to_pm_caseA":
            state = hexcols[5]
        elif datum[4].decode('utf-8') == "noMT_to_pm":
            state = hexcols[1]
        elif datum[4].decode('utf-8') == "MT_S->P_PoffMS":
            state = hexcols[2]
        elif datum[4].decode('utf-8') == "MT_S->P_PonMS":
            state = hexcols[2]
        elif datum[4].decode('utf-8') == "convergence_error":
            state = 'red'
        else:
            continue

        if datum[4][:-2].decode('utf-8') == "MT_to_pm_caseB" or datum[4][:-2].decode('utf-8') == "MT_to_pm_caseA":
            logM = float(datum[0])
            mass_ratio = float(datum[1])
            Period = float(datum[2])
            M1f = float(datum[5])
            M2f = float(datum[6])
            Pf = float(datum[7])


        logM = float(datum[0]) - logMwidth/2.0
        Period = float(datum[2]) - Pwidth/2.0

        #Add data on caseA, caseB, or convergence problems
        rectangle = patches.Rectangle((logM,Period),logMwidth,Pwidth, facecolor = state, lw = 0, ec="none")
        rect_plot = ax.add_patch(rectangle)
        if hatch != "":
            rectangle = patches.Rectangle((logM,Period),logMwidth,Pwidth, facecolor = "none", hatch = hatch, lw = 0, ec="black", color="none")
        rect_plot = ax.add_patch(rectangle)

        #if datum[4]=="ZAMS_L2_overflow":
        #    state = "black"
        #elif datum[4]=="ZAMS_RLOF":
        #    state = "0.5"
        #elif datum[4] == "off_CHE":
        #    state = hexcols[6]
        #elif datum[4] == "MT_P->S":
        #    state = hexcols[2]
        #elif datum[4] == "PISN":
        #    state = hexcols[12]
        #elif datum[4] == "MT_to_pm_caseB":
        #    state = hexcols[8]
        #elif datum[4] == "MT_to_pm_caseA":
        #    state = hexcols[5]
        #elif datum[4] == "noMT_to_pm":
        #    state = hexcols[1]
        #elif datum[4] == "MT_S->P_PoffMS":
        #    state = hexcols[2]
        #elif datum[4] == "MT_S->P_PonMS":
        #    state = hexcols[2]
        #elif datum[4] == "convergence_error":
        #    state = 'red'
        #else:
        #    continue
    nothing, = plt.plot([],[],lw=0)
    proxies = [\
            patches.Patch(color="black",lw=0),\
            patches.Patch(color="0.5",lw=0),\
            patches.Patch(color=hexcols[6], lw=0),\
            patches.Patch(color=hexcols[8], lw=0),\
            patches.Patch(color=hexcols[5], lw=0),\
            patches.Patch(color=hexcols[12], lw=0),\
            patches.Patch(color=hexcols[2], lw=0),\
            patches.Patch(color="red", lw=0),\
            nothing
            ]
    labels = [\
            "ZAMS L2OF",\
            "ZAMS RLOF",\
            "off CHE",\
            "Case B/BB",\
            "Case AB/ABB",\
            "PISN",\
            "MT before BH forms",\
            "convergence error",\
            ""
            ]

    text(1.55, 2.75, "$\\log_{10} Z="+metallicity[:-2]+"$", fontsize = 15)
    text(1.55, 2.50, "$q_\mathrm{i}="+qratio+"$", fontsize = 15)
    #text(1.4, 1.3, "$log(Z)=$"+metallicity, fontsize = 30)
    #text(1.4, 1.1, "$q_\mathrm{i}=1.0$", fontsize = 30)

    return proxies, labels

def ngal(Mbh):
    d = 2.26*(Mbh/10)**(5.0/6.0) # in Gpc
    #if Mbh>130:
    #    print 4.0/3.0*math.pi*(d*1000)**3*(2.26)**(-3)*(0.0116)/1e10, d 
    return min(1e10,4.0/3.0*math.pi*(d*1000)**3*(2.26)**(-3)*(0.0116))

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

plot_vals = [["-2.000","0.200"],["-2.500","0.200"],\
        ["-3.000","0.200"],["-3.500","0.200"],
        ["-4.000","0.200"],["-4.500","0.200"]]
plt.figure()
gs = grd.GridSpec(3,2, wspace=0.0, hspace=0.0)
for i, plot_val in enumerate(plot_vals):
    folder = plot_val[0]
    qratio = plot_val[1]
    file_name = "data_"+folder+"_"+qratio+".dat"
    ax1 = plt.subplot(gs[i])
    proxies, labels = plot_summary(ax1,0.05, 0.05, folder, qratio, file_name = file_name)

    ax1.set_ylim(0.475,3.025)
    ax1.set_xlim(1.475,2.525)
    if i == 0 or i == 2 or i == 4:
        ax1.set_ylabel("$P_\mathrm{\\rm orb, i}\;\\rm [d]$")
    if i == 4 or i == 5:
        ax1.set_xlabel("$\log_{10}M_1\;\mathrm{[M_\odot]}$")

    if i != 0 and i != 2 and i != 4:
        ax1.yaxis.set_ticklabels([])
    elif i == 0 or i== 2 or i == 4:
        ax1.yaxis.set_ticklabels(["","","1.0","1.5","2.0","2.5","3.0"])
        
    if i != 4 and i != 5:
        ax1.xaxis.set_ticklabels([])

    if i==0:
        ax1.legend(proxies, labels, loc='upper left', bbox_to_anchor=(-0.3, 1.25),
                ncol=3, scatterpoints=1,prop={'size':7})
        
plt.savefig("plots/minisummary.pdf", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                        transparent=False)
plt.clf()
plt.close(plt.gcf())
