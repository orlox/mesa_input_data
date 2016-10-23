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
from scipy.special import *
from matplotlib.colors import LogNorm

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
          'figure.subplot.bottom':0.07,
          'figure.subplot.top':0.97,
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

plt.figure()
gs = grd.GridSpec(3,1, wspace=0.05, hspace=0.3)
axes1 = plt.subplot(gs[2])
axes0 = plt.subplot(gs[0])
axes00 = plt.subplot(gs[1])

folders = ["-2.500","-3.000","-3.500","-4.000","-4.500","-5.000","-5.500","-6.000"]
#folders = ["-2.500"]
qratios = ["0.050","0.100","0.150","0.200","0.250","0.300","0.350","0.400","0.450","0.500","0.550","0.600"]
luminosities = []
spins = []
superedd = []
bh_masses = []
donor_masses = []
mass_ratios = []
cases_a = []
surfs_h1 = []
surfs_he4 = []
periods = []
weights = []
for i, folder in enumerate(folders):
    ZdivZsun_a = np.power(10,float(folder)+0.250)/0.017
    ZdivZsun_b = np.power(10,float(folder)-0.250)/0.017
    Zfrac = gammainc(0.84,ZdivZsun_a**2) - gammainc(0.84,ZdivZsun_b**2)
    print(folder, Zfrac)
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
        cases_a.append(data[5])
        periods.append(data[15])
        surfs_h1.append(data[13])
        surfs_he4.append(data[14])
        donor_masses.append(data[7])
        mass_ratios.append(data[7]/data[6])
        weights.append(data[4]*data[8]*Zfrac)

luminosities = np.power(10,np.concatenate(luminosities))*Lsun/1e39
spins = np.concatenate(spins)
cases_a = np.concatenate(cases_a)
periods = np.concatenate(periods)
surfs_h1 = np.concatenate(surfs_h1)
surfs_he4 = np.concatenate(surfs_he4)
superedd = np.concatenate(superedd)
bh_masses = np.log10(np.concatenate(bh_masses))
donor_masses = np.log10(np.concatenate(donor_masses))
mass_ratios = np.concatenate(mass_ratios)
weights = np.concatenate(weights)

sum_weights = np.sum(weights)
#weights = weights/sum_weights

print("sources per SFR:", sum_weights*0.01/3)

qbins = np.linspace(0,2)
xvals = [(qbins[i]+qbins[i+1])/2 for i in range(len(qbins)-1)]

filter_arr = np.logical_and((cases_a > 0.5), (surfs_h1 > surfs_he4))
#filter_arr = (surfs_h1 > surfs_he4)
hist, bin_edges = np.histogram(mass_ratios[filter_arr],\
        bins=qbins, weights = weights[filter_arr]*0.01/3)
for k in range(len(hist)):
    hist[k] = hist[k] / (qbins[k+1]-qbins[k])
axes00.plot(xvals, hist, color=hexcols[3], label="$q$, Case A, $X_{\\rm s} > Y_{\\rm s}$",ls="--")

filter_arr = np.logical_and((cases_a > 0.5), (surfs_h1 < surfs_he4))
#filter_arr = (surfs_h1 < surfs_he4)
hist, bin_edges = np.histogram(mass_ratios[filter_arr],\
        bins=qbins, weights = weights[filter_arr]*0.01/3)
for k in range(len(hist)):
    hist[k] = hist[k] / (qbins[k+1]-qbins[k])
axes00.plot(xvals, hist, color=hexcols[5], label="$q$, Case A, $X_{\\rm s} < Y_{\\rm s}$",ls="-.")

filter_arr = (cases_a < 0.5)
hist, bin_edges = np.histogram(mass_ratios[filter_arr],\
        bins=qbins, weights = weights[filter_arr]*0.01/3)
for k in range(len(hist)):
    hist[k] = hist[k] / (qbins[k+1]-qbins[k])
axes00.plot(xvals, hist, color=hexcols[7], label="$q$, Case AB/B",ls=":")

axes00.legend(loc=1, ncol=2, scatterpoints=1,prop={'size':7})

axes00.set_yscale("log")
axes00.set_ylim(0.001,1)
axes00.set_xlim(0,2)
axes00.set_ylabel("$dN/dq$")
axes00.set_xlabel("$q=M_2/M_{\\rm BH}$")

Mbins = 10**np.linspace(0,3,50)
xvals = [(Mbins[i]+Mbins[i+1])/2 for i in range(len(Mbins)-1)]
hist, bin_edges = np.histogram(np.power(10,bh_masses), bins=Mbins, weights = weights*0.01/3)
for k in range(len(hist)):
    hist[k] = hist[k] / (np.log10(Mbins[k+1])-np.log10(Mbins[k]))/np.log(10)
axes0.plot(xvals, hist, color=hexcols[1], label="$M_{\\rm BH}$")

filter_arr = (cases_a < 0.5)
hist, bin_edges = np.histogram(np.power(10,donor_masses[filter_arr]),\
        bins=Mbins, weights = weights[filter_arr]*0.01/3)
for k in range(len(hist)):
    hist[k] = hist[k] / (np.log10(Mbins[k+1])-np.log10(Mbins[k]))/np.log(10)
axes0.plot(xvals, hist, color=hexcols[7], label="$M_{\\rm donor}$, Case AB/B",ls=":")

filter_arr = np.logical_and((cases_a > 0.5), (surfs_h1 > surfs_he4))
#filter_arr = (surfs_h1 > surfs_he4)
hist, bin_edges = np.histogram(np.power(10,donor_masses[filter_arr]),\
        bins=Mbins, weights = weights[filter_arr]*0.01/3)
for k in range(len(hist)):
    hist[k] = hist[k] / (np.log10(Mbins[k+1])-np.log10(Mbins[k]))/np.log(10)
axes0.plot(xvals, hist, color=hexcols[3], label="$M_{\\rm donor}$, Case A, $X_{\\rm s} > Y_{\\rm s}$",ls="--")

filter_arr = np.logical_and((cases_a > 0.5), (surfs_h1 < surfs_he4))
#filter_arr = (surfs_h1 < surfs_he4)
hist, bin_edges = np.histogram(np.power(10,donor_masses[filter_arr]),\
        bins=Mbins, weights = weights[filter_arr]*0.01/3)
for k in range(len(hist)):
    hist[k] = hist[k] / (np.log10(Mbins[k+1])-np.log10(Mbins[k]))/np.log(10)
axes0.plot(xvals, hist, color=hexcols[5], label="$M_{\\rm donor}$, Case A, $X_{\\rm s} < Y_{\\rm s}$",ls="-.")

axes0.legend(loc=1, ncol=2, scatterpoints=1,prop={'size':7})

axes0.set_yscale("log")
axes0.set_xscale("log")
axes0.set_ylim(0.00005,4)
axes0.set_xlim(1,300)
axes0.set_ylabel("$M \\times dN/dM$")
axes0.set_xlabel("$M\\;[M_\\odot]$")

filter_arr = np.logical_and((cases_a > 0.5), (surfs_h1 > surfs_he4))
#filter_arr = (surfs_h1 > surfs_he4)
hist, bin_edges = np.histogram(periods[filter_arr],\
        bins=Mbins, weights = weights[filter_arr]*0.01/3)
for k in range(len(hist)):
    hist[k] = hist[k] / (np.log10(Mbins[k+1])-np.log10(Mbins[k]))/np.log(10)
axes1.plot(xvals, hist, color=hexcols[3], label="$P_{\\rm orb}$, Case A, $X_{\\rm s} > Y_{\\rm s}$",ls="--")

filter_arr = np.logical_and((cases_a > 0.5), (surfs_h1 < surfs_he4))
#filter_arr = (surfs_h1 < surfs_he4)
hist, bin_edges = np.histogram(periods[filter_arr],\
        bins=Mbins, weights = weights[filter_arr]*0.01/3)
for k in range(len(hist)):
    hist[k] = hist[k] / (np.log10(Mbins[k+1])-np.log10(Mbins[k]))/np.log(10)
left,right = bin_edges[:-1],bin_edges[1:]
X_Z = np.array([left,right]).T.flatten()
Y_Z = np.array([hist,hist]).T.flatten()
axes1.plot(xvals, hist, color=hexcols[5], label="$P_{\\rm orb}$, Case A, $X_{\\rm s} < Y_{\\rm s}$",ls="-.")

filter_arr = (cases_a < 0.5)
hist, bin_edges = np.histogram(periods[filter_arr],\
        bins=Mbins, weights = weights[filter_arr]*0.01/3)
for k in range(len(hist)):
    hist[k] = hist[k] / (np.log10(Mbins[k+1])-np.log10(Mbins[k]))/np.log(10)
axes1.plot(xvals, hist, color=hexcols[7], label="$P_{\\rm orb}$, Case AB/B",ls=":")

axes1.legend(loc=1, ncol=2, scatterpoints=1,prop={'size':8})

axes1.set_yscale("log")
axes1.set_xscale("log")
axes1.set_ylim(0.00005,4)
axes1.set_xlim(1,300)
axes1.set_ylabel("$P_{\\rm orb}\\times dN/dP_{\\rm orb}$")
axes1.set_xlabel("$P_{\\rm orb}\\;[\\rm days]$")

print("special fractions")
filter_arr = np.logical_and((donor_masses > np.log10(8)), (donor_masses < np.log10(30)))
filter_arr = np.logical_and(filter_arr, (surfs_h1 > surfs_he4))
filter_arr = np.logical_and(filter_arr, (periods<20))
filter_arr = np.logical_and(filter_arr, (cases_a>0.5))
print("donors between 8-30 msun, P<20: ", sum(weights[filter_arr])/sum(weights))
filter_arr = np.logical_and((donor_masses > np.log10(30)), (donor_masses < np.log10(70)))
filter_arr = np.logical_and(filter_arr, (cases_a>0.5))
print("donors between 30-70 msun:", sum(weights[filter_arr])/sum(weights))
filter_arr = np.logical_and((donor_masses > np.log10(30)), (donor_masses < np.log10(70)))
filter_arr = np.logical_and(filter_arr, (surfs_h1 > surfs_he4))
filter_arr = np.logical_and(filter_arr, (cases_a>0.5))
print("donors between 30-70 msun (H-rich): ", sum(weights[filter_arr])/sum(weights))
filter_arr = np.logical_and((donor_masses > np.log10(30)), (donor_masses < np.log10(70)))
filter_arr = np.logical_and(filter_arr, (surfs_h1 < surfs_he4))
filter_arr = np.logical_and(filter_arr, (cases_a>0.5))
print("donors between 30-70 msun (H-poor): ", sum(weights[filter_arr])/sum(weights))
filter_arr = (cases_a < 0.5)
print("case B:", sum(weights[filter_arr])/sum(weights))
filter_arr = (mass_ratios < 0.5)
print("q<0.5:", sum(weights[filter_arr])/sum(weights))

plt.savefig("plots/mass_dist.pdf", dpi=300, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                        transparent=False)
