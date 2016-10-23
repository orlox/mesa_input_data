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
          'figure.subplot.right':0.9}


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

range_integral = 0.0194213341446518*(2.56+0.301)
limited_range = 0.004105169206845768*(0.5+0.301)

folders = ["-2.500","-3.000","-3.500","-4.000","-4.500","-5.000","-5.500","-6.000"]
Zvals = [-2.50,-3.00,-3.50,-4.00,-4.50,-5.00,-5.50,-6.00]
BHNS = np.zeros(len(folders))
BHNS_disrupted = np.zeros(len(folders))
BHNS_merger = np.zeros(len(folders))
BB = np.zeros(len(folders))
BHBH = np.zeros(len(folders))
BHBH_disrupted = np.zeros(len(folders))
BHBH_merger = np.zeros(len(folders))
PISN = np.zeros(len(folders))
for i, folder in enumerate(folders):
    file_name = "kick_"+folder+".dat"
    data = np.genfromtxt(file_name, dtype=None, names=True, skip_header=1)
    print file_name
    for datum in data:
        if float(datum[2])<0.49:
            continue
        logM = float(datum[0])
        dlm = 0.05
        Period = float(datum[2])
        dlp = np.log10(Period+0.025) - np.log10(Period-0.025)
        #dlp = ((Period+0.025)*np.log(Period+0.025)-(Period+0.025))/np.log(10) - \
        #        ((Period-0.025)*np.log(Period-0.025)-(Period-0.025))/np.log(10)
        dq = 0.05
        extra_weight = np.power(np.power(10,logM),-1.35)*dlm*dlp*dq
        if datum[6] > 0:
            BB[i] += extra_weight
        elif datum[4] < 12:
            BHNS[i] += (1-datum[8]-datum[7])*extra_weight
            BHNS_disrupted[i] += datum[8]*extra_weight
            BHNS_merger[i] += datum[7]*extra_weight
        elif datum[4] < 60:
            BHBH[i] += (1-datum[10]-datum[9])*extra_weight
            BHBH_disrupted[i] += datum[10]*extra_weight
            BHBH_merger[i] += datum[9]*extra_weight
        else:
            PISN[i] += extra_weight
    sum_remnants = (BB[i]\
            +BHNS[i]+BHNS_disrupted[i]+BHNS_merger[i]\
            +BHBH[i]+BHBH_disrupted[i]+BHBH_merger[i]\
            +PISN[i])
    BHNS[i] /= sum_remnants
    BHNS_disrupted[i] /= sum_remnants
    BHNS_merger[i] /= sum_remnants
    BB[i] /= sum_remnants
    BHBH[i] /= sum_remnants
    BHBH_disrupted[i] /= sum_remnants
    BHBH_merger[i] /= sum_remnants
    PISN[i] /= sum_remnants

print "BB,BHNS,BHNSd,BHNSm,BHBH,BHBHd,BHBHm,PISN"
for i in range(len(BB)):
    print BB[i], BHNS[i],BHNS_disrupted[i],BHNS_merger[i],BHBH[i],BHBH_disrupted[i],BHBH_merger[i],PISN[i]

plt.plot(Zvals,np.log10(BB+1e-10),color="0.5",label="Case MBB")
plt.plot(Zvals,np.log10(BHNS+1e-10),color=hexcols[1],label="wide BH-NS")
plt.plot(Zvals,np.log10(BHNS_merger+1e-10), color=hexcols[1], ls="--",label="merging BH-NS")
plt.plot(Zvals,np.log10(BHNS_disrupted+1e-10), color=hexcols[1], ls =":",label="disrupted BH-NS")
plt.plot(Zvals,np.log10(BHBH+1e-10), color=hexcols[3],label="wide BH-BH")
plt.plot(Zvals,np.log10(BHBH_merger+1e-10), color=hexcols[3],ls="--",label="merging BH-BH")
plt.plot(Zvals,np.log10(BHBH_disrupted+1e-10), color=hexcols[3],ls=":",label="disrupted BH-BH")
plt.plot(Zvals,np.log10(PISN+1e-10), color = hexcols[5],label="PISN")

plt.legend(loc=3)
plt.gca().set_ylim([-3.5,0])

plt.savefig("post_kick.pdf", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                        transparent=False)


#    import colormaps as cmaps
#    plt.register_cmap(name='viridis', cmap=cmaps.viridis)
#    plt.set_cmap(cmaps.viridis)
#
#    #gs = grd.GridSpec(2, 2, height_ratios=[1,2], width_ratios=[20,1], wspace=0.1, hspace=0.1)
#    gs = grd.GridSpec(1, 1)
#
#    axarr = [plt.subplot(gs[0]),plt.subplot(gs[2]),plt.subplot(gs[3])]
#
#    axarr[0].set_ylabel("relative number")
#    axarr[0].xaxis.set_major_locator(plt.NullLocator())
#    n, bins, patches = axarr[0].hist( [bh_masses_A, bh_masses_B], weights=[bh_weights_A,bh_weights_B],\
#            histtype='barstacked', stacked = True, color =[hexcols[1],hexcols[4]], label=["Case A","Case B"],\
#            rwidth = 1, bins=np.arange(1.4,2.7,0.05),lw=0)
#    legend = axarr[0].legend(loc="best",title="$\log\;Z="+folder+"$\n"+"{0:.3}".format(sum_weights*0.01)+" per unit SFR")
#    axarr[0].set_xlim([1.4,2.6])
#    axarr[0].set_ylim([0,0.4])
#    plt.setp(legend.get_title(),fontsize='7')
#
#    axarr[1].set_ylabel("$\log$ age [yr]")
#    axarr[1].set_xlabel("$\log\;M_{\mathrm{BH}}\;\mathrm{[M_\odot]}$")
#    H, xedges, yedges, img = axarr[1].hist2d( hist_data_bhmass, hist_data_times, weights=hist_data_weights,\
#            bins=[np.arange(1.4,2.7,0.05),np.arange(6.3,8.15,0.1)], normed = True)
#    axarr[1].set_xlim([1.4,2.6])
#
#    norm = mpl.colors.Normalize(vmin=0, vmax=1)
#    mpl.colorbar.ColorbarBase(axarr[2], norm=norm, orientation='vertical')
#
#    plt.savefig("plots/mt_times_"+folder+".pdf", dpi=None, facecolor='w', edgecolor='w',
#                    orientation='portrait', papertype=None, format=None,
#                            transparent=False)
