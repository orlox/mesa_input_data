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
from mesa import *
import mmap

params = {'backend': 'pdf',
          'figure.figsize': [4.3, 6.0],
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
          'figure.subplot.bottom':0.1,
          'figure.subplot.top':0.9,
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

def generate_cols(history_file):
    """
    Returns a dictionary with column numbers (starting from zero) for the mesa
    output file history_file
    """
    cols = {}
    file = open(history_file, "r")
    for i in range(4):
        file.readline()
    nums = file.readline().split()
    names = file.readline().split()
    for i, name in enumerate(names):
        cols[name] = int(nums[i])-1

def plot_summary(ax,masses_donor, periods, qratio, folder ='.'):
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
    folder = folder + "/"

    #solve boundaries in period for the plot
    boundaries_period = zeros(len(periods)+1)
    for i in range(1,len(boundaries_period)-1):
        boundaries_period[i] = (float(periods[i-1]) + float(periods[i])) / 2
    boundaries_period[0] = 2*float(periods[0]) - boundaries_period[1]
    boundaries_period[-1] = 2*float(periods[-1]) - boundaries_period[-2]
    #solve boundaries in accretor mass for the plot
    boundaries_mass = zeros(len(masses_donor)+1)
    for i in range(1,len(boundaries_mass)-1):
        boundaries_mass[i] = (float(masses_donor[i-1]) + float(masses_donor[i])) / 2
    boundaries_mass[0] = 2*float(masses_donor[0]) - boundaries_mass[1]
    boundaries_mass[-1] = 2*float(masses_donor[-1]) - boundaries_mass[-2]

    masses_i = np.zeros(len(masses_donor))
    periods_i = np.zeros(len(periods))

    for i, period in enumerate(periods):
        periods_i[i] = float(period)
    for i, mass in enumerate(masses_donor):
        masses_i[i] = float(mass)

    lperiod_f = np.ones((len(masses_donor),len(periods)))*-1e50
    lmass_f = np.ones((len(masses_donor),len(periods)))*-1e50

    masses_end = []
    masses_end2 = []
    periods_end = []
    pointslM = []
    pointslP = []
    contacts = []
    kerr_params = []

    for i, period in enumerate(periods):
        for j, mass_donor in enumerate(masses_donor):
            #try:
                state = ""
                hatch = ""
                folder_name = folder+mass_donor+"_"+qratio+"_"+period+"/"
                if not os.path.isdir(folder_name):
                   print "folder does not exist", folder_name
                   continue
                outdata = os.popen('tail -n 40 '+folder_name+'out.txt').read()
                outfile = open(folder_name+'out.txt')
                s = mmap.mmap(outfile.fileno(), 0, access=mmap.ACCESS_READ)
                if outdata.find('DATE:') == -1:
                    continue
                if s.rfind('model is overflowing at ZAMS') != -1:
                    hatch = "xxxx"
                if outdata.find('Terminate due to overflow of L2 at ZAMS') != -1:
                    state = "black"
                elif float(period)<0.9 and float(mass_donor)>2.3:
                    state = "black"
                elif outdata.find('Terminate due to secondary not evolving homogeneously') != -1:
                    state = hexcols[6]
                elif outdata.find('Terminate due to primary not evolving homogeneously') != -1:
                    state = hexcols[6]
                elif outdata.find('termination code: Terminate because of L2 overflow') != -1:
                    state = hexcols[2]
                elif outdata.find('Terminate due to helium depletion') != -1:
                    state = hexcols[1]
                    binary_history = history_data(folder_name, slname = "binary_history.data", clean_starlog=False)
                    star_history = history_data(folder_name+"/LOGS1", slname = "history.data", clean_starlog=False)
                    masses_end.append(binary_history.get("star_1_mass")[-1])
                    masses_end2.append(binary_history.get("star_2_mass")[-1])
                    periods_end.append(binary_history.get("period_days")[-1])
                    lperiod_f[j,i] = np.log10(binary_history.get("period_days")[-1])
                    lmass_f[j,i] = np.log10(2*binary_history.get("star_1_mass")[-1])
                    kerr_params.append(\
                            np.power(10,star_history.get("log_total_angular_momentum")[-1])\
                            *clight/((binary_history.get("star_1_mass")[-1]*msun)**2*cgrav)
                            )
                    if kerr_params[-1]<0.2 and np.log10(2*masses_end[-1])<2.4:
                        print "LOW KERR PARAMETER!!!", kerr_params[-1], np.log10(2*masses_end[-1])
                    if kerr_params[-1]<0.4 and np.log10(2*masses_end[-1])<2.2:
                        print "LOW KERR PARAMETERB!!!"
                    pointslM.append(float(mass_donor))
                    pointslP.append(np.log10(float(period)))
                    contact = 0
                    if hatch == "xxxx":
                        contact = 2
                    else:
                        overflow_1 = binary_history.get("rl_relative_overflow_1")
                        overflow_2 = binary_history.get("rl_relative_overflow_2")
                        for k in range(len(overflow_1)):
                            if overflow_1[k]>0 and overflow_2[k]>0:
                                hatch = "////"
                                contact = 1
                                break
                    contacts.append(contact)
                else:
                    state = "red"

                print folder_name

                #Add data on caseA, caseB, or convergence problems
                rectangle = patches.Rectangle((boundaries_mass[j],boundaries_period[i]),\
                      boundaries_mass[j+1]-boundaries_mass[j],\
                      boundaries_period[i+1]-boundaries_period[i], facecolor = state, lw = 0, ec="none")
                rect_plot = ax.add_patch(rectangle)
                if hatch != "":
                    rectangle = patches.Rectangle((boundaries_mass[j],boundaries_period[i]),\
                          boundaries_mass[j+1]-boundaries_mass[j],\
                          boundaries_period[i+1]-boundaries_period[i], facecolor = "none", hatch = hatch, lw = 0, ec="black", color="none")
                rect_plot = ax.add_patch(rectangle)

    proxies = [\
            patches.Patch(color="black",lw=0),\
            patches.Patch(color=hexcols[6], lw=0),\
            patches.Patch(facecolor="none",lw=0, hatch ="xxxx"),\
            patches.Patch(color=hexcols[2], lw=0),\
            patches.Patch(facecolor="none",lw=0, hatch ="////"),\
            patches.Patch(color=hexcols[1], lw=0),\
            ]
    labels = [\
            "ZAMS L2OF",\
            "Off CHE",\
            "ZAMS RLOF",\
            "L2 overflow",\
            "contact MS",\
            "Double he star",\
            ]

    text(1.4, 2.75, "$Z=Z_\\odot/50$", fontsize = 30)
    text(1.4, 2.45, "$q_\mathrm{i}=1.0$", fontsize = 30)
    
    plt.ylim(boundaries_period[0],boundaries_period[-1])
    plt.xlim(boundaries_mass[0],boundaries_mass[-1])

    return [boundaries_mass[0],boundaries_mass[-1]], [boundaries_period[0],boundaries_period[-1]], \
            proxies, labels, masses_end, masses_end2, periods_end, pointslM, pointslP, contacts, kerr_params

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

lg_masses_donor =    ["1.4","1.5","1.6","1.7","1.8","1.9","2.0","2.1","2.2","2.3","2.4","2.5","2.6","2.7"]
periods = ["0.50","0.60","0.70","0.80","0.90","1.00","1.10","1.20"\
        ,"1.30","1.40","1.50","1.60","1.70","1.80","1.90"\
        ,"2.00","2.10","2.20","2.30","2.40","2.50","2.60","2.70","2.80","2.90","3.00"]
#lg_masses_donor =    ["2.2","2.4"]
#periods = ["1.50","1.70"]

qratios = ["1.0","0.9","0.8"]
#qratios = ["1.0"]

folders = ["Z50","Z20","Z10","Z04"]
folders = ["Z50","Z20","Z10"]
#folders = ["Z50"]

masses_Z = {}
masses2_Z = {}
periods_Z = {}
pointslM_Z = {}
pointslP_Z = {}
contacts_Z = {}
kerr_params_Z = {}
for folder in folders:
    for qratio in qratios:
        if qratio != "1.0" and folder != "Z50":
            continue
        folder_name = "double_bh_"+ folder
        plt.figure()
        f, (ax1) = plt.subplots(1, 1, sharex=True)
        xlims, ylims, proxies, labels, mass_Z, mass2_Z, period_Z, pointlM_Z, pointlP_Z, contact_Z, kerr_param_Z = \
                plot_summary(ax1,lg_masses_donor, periods, qratio, folder = folder_name)
        if qratio == "1.0":
            masses_Z[folder] = mass_Z
            periods_Z[folder] = period_Z
            pointslM_Z[folder] = pointlM_Z
            pointslP_Z[folder] = pointlP_Z
            contacts_Z[folder] = contact_Z
            kerr_params_Z[folder] = kerr_param_Z
            print kerr_param_Z
        else:
            masses_Z[folder+qratio] = mass_Z
            masses2_Z[folder+qratio] = mass2_Z
            periods_Z[folder+qratio] = period_Z
            pointslM_Z[folder+qratio] = pointlM_Z
            pointslP_Z[folder+qratio] = pointlP_Z
            contacts_Z[folder+qratio] = contact_Z
        
        ax1.set_ylim(ylims[0],ylims[1])
        ax1.set_ylabel("$\mathrm{P_\mathrm{i}\;[d]}$")
        ax1.set_xlabel("$\mathrm{\log\;M_1\;\mathrm{[M_\odot]}}$")
        
        ax1.legend(proxies, labels, loc='upper center', bbox_to_anchor=(0.5, 1.1),
                ncol=3, scatterpoints=1,prop={'size':7})
        
        plt.savefig("summary_"+qratio+"_"+folder+".pdf", dpi=None, facecolor='w', edgecolor='w',
                        orientation='portrait', papertype=None, format=None,
                                transparent=False)

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
          'figure.subplot.right':1.0}

mpl.rcParams.update(params)

plt.figure()

xpoints = {}
ypoints = {}
tbin = 0.005
tvals = np.arange(0,1,tbin)
hist_data = {\
        "Z50":np.zeros(len(tvals)),\
        "Z20":np.zeros(len(tvals)),\
        "Z10":np.zeros(len(tvals)),\
        "Z04":np.zeros(len(tvals))\
        }

lims = {\
        "Z50": np.array([[1.81,0.1],[2.75,1.25]]),\
        "Z20": np.array([[1.75,0.25],[2.45,1.5]]),\
        "Z10": np.array([[1.65,0.4],[2.2,1.65]]),\
        "Z04": np.array([[1.72,1],[1.97,2.5]])
        }
for j, folder in enumerate(folders):
    xpoints[folder] = (1.0-tvals)*lims[folder][0][0]+tvals*lims[folder][1][0]
    ypoints[folder] = (1.0-tvals)*lims[folder][0][1]+tvals*lims[folder][1][1]

dlm = 0.00025
dlp = 0.00025
#dlm = 0.001
#dlp = 0.001
#dlm = 0.1
#dlp = 0.1
#lmcoords = np.arange(1.4,2.8,dm)
#lpcoords = np.arange(-0.301, 0.5, dlp)
lmcoords, lpcoords = np.mgrid[1.4:2.9:dlm,-0.301:0.5:dlp]

dlm_chirp = 0.01
chirp_vals = np.arange(0.0,3.0,dlm_chirp)
chirp_data = np.zeros((len(folders),len(chirp_vals),3))
chirp_LIGO_data = np.zeros((len(folders),len(chirp_vals),3))
chirp_LIGO_data_low = np.zeros((len(folders),len(chirp_vals),3))
chirp_LIGO_data_high = np.zeros((len(folders),len(chirp_vals),3))

#dt_merger = 0.1
#merger_vals = np.arange(0.0,15,dt_merger)
#merger_data = np.zeros((len(folders),len(merger_vals),3))
dt_merger = 0.01
merger_vals = np.arange(-1,1.5,dt_merger)
merger_data = np.zeros((len(folders),len(merger_vals),3))

Z_weights = [0.035,0.04,0.1,1.0]
Z_weights_low = [0.0004,0.0025,0.01,1.0]
Z_weights_high = [0.068,0.052,0.086,1.0]
weights = np.power(np.power(10,lmcoords[:,0]),-1.35)*dlm*dlp
print "test weights", sum( np.power(np.power(10,lmcoords[:,0]),-1.35)*dlm), 0.004105169206845768
range_integral = 0.0194213341446518*(2.56+0.301)
limited_range = 0.004105169206845768*(0.5+0.301)
rates = np.zeros((len(folders),5))
rates_LIGO = np.zeros((len(folders),5))
sum_check = 0
tmerge_min = 1e10
for j, folder in enumerate(folders):
    print "getting mass data for",folder
    #mass_val = griddata(pointslM_Z[folder], pointslP_Z[folder], np.log10(2*np.array(masses_Z[folder])), lmcoords,lpcoords, interp='nn')
    mass_vals = griddata(zip(pointslM_Z[folder], pointslP_Z[folder]), np.log10(2*np.array(masses_Z[folder])), (lmcoords,lpcoords), method='linear')
    print "getting period data for",folder
    period_vals = griddata(zip(pointslM_Z[folder], pointslP_Z[folder]), np.log10(np.array(periods_Z[folder])), (lmcoords,lpcoords), method='linear')
    print "fraction of interpolated progenitors:", np.count_nonzero(~np.isnan(mass_vals)), "/", mass_vals.size, "=", \
            float(np.count_nonzero(~np.isnan(mass_vals)))/float(mass_vals.size)
    contact_vals = griddata(zip(pointslM_Z[folder], pointslP_Z[folder]), contacts_Z[folder], (lmcoords,lpcoords), method='linear')
    for i in range(len(lmcoords)):
        for k in range(len(lpcoords[0])):
            sum_check += weights[i]
            if math.isnan(mass_vals[i,k]):
                continue
            mass_val = mass_vals[i,k]
            if math.isnan(period_vals[i,k]):
                print "something is wrong!"
                continue
            P_val = period_vals[i,k]

            merges = False
            if P_val < 3.0/8.0*(np.log10(13.8)-2.276+5.0/3.0*mass_val):
                merges = True

            merge_time = pow(10, 8.0/3.0*P_val+2.276-5.0/3.0*mass_val)
            tmerge_min = min(tmerge_min, merge_time)

            Mbh = np.power(10,mass_val)/2.0

            mass_range = 1
            if mass_val < 2.079:
                mass_range = 0
            elif mass_val > 2.415:
                mass_range = 2
            if mass_range == 1:
                rates[j,2] += weights[i]
                rates_LIGO[j,2] += weights[i]*0.01*(1.0/3.0)*0.2*ngal(Mbh)
            elif mass_range == 0:
                if merges:
                    rates[j,0] += weights[i]
                    rates_LIGO[j,0] += weights[i]*0.01*(1.0/3.0)*0.2*ngal(Mbh)
                else:
                    rates[j,1] += weights[i]
                    rates_LIGO[j,1] += weights[i]*0.01*(1.0/3.0)*0.2*ngal(Mbh)
            elif mass_range == 2:
                if merges:
                    rates[j,3] += weights[i]
                    rates_LIGO[j,3] += weights[i]*0.01*(1.0/3.0)*0.2*ngal(Mbh)
                else:
                    rates[j,4] += weights[i]
                    rates_LIGO[j,4] += weights[i]*0.01*(1.0/3.0)*0.2*ngal(Mbh)

            if (mass_range == 0 or mass_range == 2) and merges:
                lMchirp = np.log10(0.87*Mbh)
                t = int(lMchirp/dlm_chirp)
                if t< 0:
                    print "ERROR!!!"
                elif t>len(chirp_vals)-1:
                    print "ERROR!!!"
                had_contact = 0
                if contact_vals[i,k] > 0.5 and contact_vals[i,k] < 1.5:
                    had_contact = 1
                elif contact_vals[i,k] >= 1.5:
                    had_contact = 2
                chirp_data[j,t,had_contact] += weights[i]*Z_weights[j]
                chirp_LIGO_data[j,t,had_contact] += weights[i]*0.01*(1.0/3.0)*0.2*ngal(Mbh)*Z_weights[j]
                chirp_LIGO_data_low[j,t,had_contact] += weights[i]*0.01*(1.0/3.0)*0.2*ngal(Mbh)*Z_weights_low[j]
                chirp_LIGO_data_high[j,t,had_contact] += weights[i]*0.01*(1.0/3.0)*0.2*ngal(Mbh)*Z_weights_high[j]

            if (mass_range == 0 or mass_range == 2):
                #t = int(merge_time/dt_merger)
                #if t< 0:
                #    print "ERROR!!!"
                t = int((np.log10(merge_time)+1)/dt_merger)
                if t< 0:
                    print "ERROR!!!"
                elif t>len(merger_vals)-1:
                    continue
                had_contact = 0
                if contact_vals[i,k] > 0.5 and contact_vals[i,k] < 1.5:
                    had_contact = 1
                elif contact_vals[i,k] >= 1.5:
                    had_contact = 2
                merger_data[j,t,had_contact] += weights[i]*Z_weights[j]



            vec_x0 = [mass_val,P_val]
            t = np.dot(lims[folder][0], lims[folder][1])\
                    +np.dot(lims[folder][0], vec_x0)\
                    -np.dot(lims[folder][0], lims[folder][0])\
                    -np.dot(lims[folder][1], vec_x0)
            t = t/(2*np.dot(lims[folder][0], lims[folder][1])\
                    -np.dot(lims[folder][0], lims[folder][0])\
                    -np.dot(lims[folder][1], lims[folder][1]))
            t = int(t/tbin)
            if t< 0:
                t = 0
            elif t>len(tvals)-1:
                t = len(tvals)-1
            hist_data[folder][t] += weights[i]

    print "checking integration weights", sum_check, limited_range, sum_check / limited_range
    print "tmerge_min is", tmerge_min
        
#for j, folder in enumerate(folders):
#    total = sum(systems)

print rates/range_integral
print rates/range_integral*0.01*(1.0/3.0)*0.2*1e6
print rates_LIGO/range_integral

# histogram the data
xyrange = [[1.35,2.85],[-0.2,1.9]]
#xyrange = [[1.35,2.85],[-0.2,3.5]]
#hh, locx, locy = scipy.histogram2d(xdat, ydat, range=xyrange, bins=[40,25])
#plt.imshow(np.flipud(hh.T),cmap='hot_r',extent=np.array(xyrange).flatten(), interpolation='none',aspect='auto')

plt.gca().add_patch(
    patches.Rectangle(
        (2.07, -1.0),   # (x,y)
        2.415-2.07,          # width
        5,          # height
    color = colorscale(hexcols[6],2), alpha = 0.7, zorder = -1000)
)

text(2.075, -0.05, "PISN", fontsize = 22, color = hexcols[6])

#plt.scatter(xpoints["Z04"], ypoints["Z04"], c=hist_data["Z04"]/max(hist_data["Z04"]), s=250, marker = "o", linewidth = 0, edgecolor = "none", cmap ="hot_r", zorder = -50)
plt.scatter(xpoints["Z10"], ypoints["Z10"], c=hist_data["Z10"]/max(hist_data["Z10"]), s=250, marker = "o", linewidth = 0, edgecolor = "none", cmap ="hot_r", zorder = -50)
plt.scatter(xpoints["Z20"], ypoints["Z20"], c=hist_data["Z20"]/max(hist_data["Z20"]), s=250, marker = "o", linewidth = 0, edgecolor = "none", cmap ="hot_r", zorder = -50)
plt.scatter(xpoints["Z50"], ypoints["Z50"], c=hist_data["Z50"]/max(hist_data["Z50"]), s=250, marker = "o", linewidth = 0, edgecolor = "none", cmap ="hot_r", zorder = -50)

cb = plt.colorbar(ticks=[0.0, 0.5, 1.0])

#plt.scatter(np.log10(2*np.array(masses_Z["Z04"])),np.log10(periods_Z["Z04"]), s=30, marker = "v", facecolor = hexcols[7], edgecolor =colorscale(hexcols[5],0.5), linewidth = 1, label = "$Z_\\odot/4$")
plt.scatter(np.log10(2*np.array(masses_Z["Z10"])),np.log10(periods_Z["Z10"]), s=30, marker = "^", facecolor = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), linewidth = 1, label = "$Z_\\odot/10$")
plt.scatter(np.log10(2*np.array(masses_Z["Z20"])),np.log10(periods_Z["Z20"]), s=30, marker = "o", facecolor = hexcols[5], edgecolor =colorscale(hexcols[7],0.5), linewidth = 1, label = "$Z_\\odot/20$")
plt.scatter(np.log10(2*np.array(masses_Z["Z50"])),np.log10(periods_Z["Z50"]), s=30, marker = "s", facecolor = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), linewidth = 1, label = "$Z_\\odot/50$")

mrange = np.array([1.7,10])

plt.plot(mrange, 3.0/8.0*(np.log10(1)-2.276+5.0/3.0*mrange),'--',color="0.5", zorder = -10)
plt.plot(mrange, 3.0/8.0*(np.log10(4)-2.276+5.0/3.0*mrange),'--',color="0.5", zorder = -10)
plt.plot(mrange, 3.0/8.0*(np.log10(13.8)-2.276+5.0/3.0*mrange),'--',color="0.5", zorder = -10)

text(1.38, 0.57, "$13.8\;\mathrm{Gyr}$", fontsize = 10, rotation = 16, color = "0.5")
text(1.44, 0.35, "$4\;\mathrm{Gyr}$", fontsize = 10, rotation = 16, color = "0.5")
text(1.47, 0.13, "$1\;\mathrm{Gyr}$", fontsize = 10, rotation = 16, color = "0.5")

#prange = np.array([-10,10])
#
#plt.plot(3.0/5.0*(-np.log10(1)+2.276+8.0/3.0*prange), prange,'--',color="0.5",zorder = -1000)
#plt.plot(3.0/5.0*(-np.log10(4)+2.276+8.0/3.0*prange), prange,'--',color="0.5",zorder = -1000)
#plt.plot(3.0/5.0*(-np.log10(13.8)+2.276+8.0/3.0*prange), prange,'--',color="0.5", zorder = -1000)

plt.legend(loc=4, scatterpoints = 1)
plt.gca().set_xlabel("$\log(M_1+M_2)\;\mathrm{[M_\odot]}$")
plt.gca().set_ylabel("$\log\;P_{\mathrm{f}}\;\mathrm{[d]}$")
plt.gca().set_xlim(xyrange[0][0],xyrange[0][1])
plt.gca().set_ylim(xyrange[1][0],xyrange[1][1])

plt.savefig("mass_period.pdf", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                        transparent=False)

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
          'figure.subplot.right':0.95}

mpl.rcParams.update(params)

plt.figure()

plt.gca().add_patch(
    patches.Rectangle(
        (2.07, -1.0),   # (x,y)
        2.415-2.07,          # width
        5,          # height
    color = colorscale(hexcols[6],2), alpha = 0.7, zorder = -1000)
)

text(2.11, 0.87, "PISN", fontsize = 22, color = hexcols[6])

#plt.scatter(np.log10(2*np.array(masses_Z["Z04"])),np.log10(periods_Z["Z04"]), s=30, marker = "v", facecolor = hexcols[7], edgecolor =colorscale(hexcols[5],0.5), linewidth = 1, label = "$Z_\\odot/4$")

Zlabels = ["$Z_\\odot/10$","$Z_\\odot/20$","$Z_\\odot/50$"]
Zcolors = [hexcols[3], hexcols[5], hexcols[1]]
Zstrs = ["Z10","Z20","Z50"]
Zmarkers = ["^","o","s"]

for i, Zstr in enumerate(Zstrs):
    Mmasked = np.ma.masked_where(np.log10(periods_Z[Zstr]) > 3.0/8.0*(np.log10(13.8)-2.276+5.0/3.0*np.log10(2*np.array(masses_Z[Zstr])))\
            , masses_Z[Zstr])
    Mmasked = np.ma.compressed(Mmasked)
    KERRmasked = np.ma.masked_where(np.log10(periods_Z[Zstr]) > 3.0/8.0*(np.log10(13.8)-2.276+5.0/3.0*np.log10(2*np.array(masses_Z[Zstr])))\
            , kerr_params_Z[Zstr])
    KERRmasked = np.ma.compressed(KERRmasked)
    plt.scatter(np.log10(2*np.array(Mmasked)),np.minimum(1.0,KERRmasked),\
            s=30, marker = Zmarkers[i], facecolor = Zcolors[i], edgecolor =colorscale(Zcolors[i],0.5), linewidth = 1, label = Zlabels[i])
    Mmasked = np.ma.masked_where(np.log10(periods_Z[Zstr]) < 3.0/8.0*(np.log10(13.8)-2.276+5.0/3.0*np.log10(2*np.array(masses_Z[Zstr])))\
            , masses_Z[Zstr])
    Mmasked = np.ma.compressed(Mmasked)
    KERRmasked = np.ma.masked_where(np.log10(periods_Z[Zstr]) < 3.0/8.0*(np.log10(13.8)-2.276+5.0/3.0*np.log10(2*np.array(masses_Z[Zstr])))\
            , kerr_params_Z[Zstr])
    KERRmasked = np.ma.compressed(KERRmasked)
    plt.scatter(np.log10(2*np.array(Mmasked)),np.minimum(1.0,KERRmasked),\
            s=30, marker = Zmarkers[i], facecolor = Zcolors[i], edgecolor ="r", linewidth = 1)
    #Mmasked = np.ma.masked_where(np.log10(periods_Z[Zstr]) > 3.0/8.0*(np.log10(13.8)-2.276+5.0/3.0*np.log10(2*np.array(masses_Z[Zstr])))\
    #        or (np.array(masses_Z[Zstr]) > 60 and np.array(masses_Z[Zstr]) < 130), masses_Z[Zstr])
    #Mmasked = np.ma.compressed(Mmasked)
    #KERRmasked = np.ma.masked_where(np.log10(periods_Z[Zstr]) > 3.0/8.0*(np.log10(13.8)-2.276+5.0/3.0*np.log10(2*np.array(masses_Z[Zstr])))\
    #        or (np.array(masses_Z[Zstr]) > 60 and np.array(masses_Z[Zstr]) < 130), kerr_params_Z[Zstr])
    #KERRmasked = np.ma.compressed(KERRmasked)
    #plt.scatter(np.log10(2*np.array(Mmasked)),np.minimum(1.0,KERRmasked),\
    #        s=30, marker = Zmarkers[i], facecolor = Zcolors[i], edgecolor =colorscale(Zcolors[i],0.5), linewidth = 1, label = Zlabels[i])
    #Mmasked = np.ma.masked_where(np.log10(periods_Z[Zstr]) < 3.0/8.0*(np.log10(13.8)-2.276+5.0/3.0*np.log10(2*np.array(masses_Z[Zstr])))\
    #        and (np.array(masses_Z[Zstr]) < 60 or np.array(masses_Z[Zstr]) > 130), masses_Z[Zstr])
    #Mmasked = np.ma.compressed(Mmasked)
    #KERRmasked = np.ma.masked_where(np.log10(periods_Z[Zstr]) < 3.0/8.0*(np.log10(13.8)-2.276+5.0/3.0*np.log10(2*np.array(masses_Z[Zstr])))\
    #        and (np.array(masses_Z[Zstr]) < 60 or np.array(masses_Z[Zstr]) > 130), kerr_params_Z[Zstr])
    #KERRmasked = np.ma.compressed(KERRmasked)
    #plt.scatter(np.log10(2*np.array(Mmasked)),np.minimum(1.0,KERRmasked),\
    #        s=30, marker = Zmarkers[i], facecolor = Zcolors[i], edgecolor ="r", linewidth = 1)

#prange = np.array([-10,10])
#
#plt.plot(3.0/5.0*(-np.log10(1)+2.276+8.0/3.0*prange), prange,'--',color="0.5",zorder = -1000)
#plt.plot(3.0/5.0*(-np.log10(4)+2.276+8.0/3.0*prange), prange,'--',color="0.5",zorder = -1000)
#plt.plot(3.0/5.0*(-np.log10(13.8)+2.276+8.0/3.0*prange), prange,'--',color="0.5", zorder = -1000)

plt.legend(loc=1, scatterpoints = 1)
plt.gca().set_xlabel("$\log(M_1+M_2)\;\mathrm{[M_\odot]}$")
plt.gca().set_ylabel("$min(1,Jc/M^2G)$")
plt.gca().set_xlim(xyrange[0][0],xyrange[0][1])
plt.gca().set_ylim(0,1.05)

plt.savefig("mass_kerr.pdf", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                        transparent=False)

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
          'figure.subplot.right':0.95}

mpl.rcParams.update(params)

plt.figure()

plt.scatter(np.log10(np.array(masses_Z["Z500.9"])+np.array(masses2_Z["Z500.9"])),np.array(masses2_Z["Z500.9"])/np.array(masses_Z["Z500.9"]), \
        s=40, marker = "s", facecolor = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), linewidth = 1, label = "$q_\mathrm{i}=0.9$")
plt.scatter(np.log10(np.array(masses_Z["Z500.8"])+np.array(masses2_Z["Z500.8"])),np.array(masses2_Z["Z500.8"])/np.array(masses_Z["Z500.8"]), \
        s=40, marker = "^", facecolor = hexcols[7], edgecolor =colorscale(hexcols[7],0.5), linewidth = 1, label = "$q_\mathrm{i}=0.8$")

plt.gca().add_patch(
    patches.Rectangle(
        (2.07, -1.0),   # (x,y)
        2.41-2.07,          # width
        5,          # height
    color = colorscale(hexcols[6],2), alpha = 0.7, zorder = -1000)
)

text(2.2, 0.88, "PISN", fontsize = 25, color = hexcols[6])

legend = plt.legend(loc=3, scatterpoints = 1, title = "$Z=Z_\\odot/50$")
legend.get_title().set_fontsize('14')
#plt.gca().set_xlim(xyrange[0][0],xyrange[0][1])
plt.gca().set_xlabel("$\log\;(M_1+M_2)$")
plt.gca().set_ylabel("$q=M_2/M_1$")
#plt.gca().set_xlim(1.35,2.85)
#plt.gca().set_ylim(-0.2,1.9)

plt.savefig("mass_q.pdf", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                        transparent=False)

proxies = [\
        patches.Patch(facecolor="none",lw=0, hatch ="////"),\
        patches.Patch(facecolor="none",lw=0, hatch ="xxxx"),\
        patches.Patch(color=hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2),\
        patches.Patch(color=hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2),\
        patches.Patch(color=hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2),\
        ]
labels = [\
        "Contact during the MS",\
        "Contact at ZAMS",\
        "$Z_\\odot/10$",\
        "$Z_\\odot/20$",\
        "$Z_\\odot/50$",\
        ]


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
          'figure.subplot.right':0.95}

mpl.rcParams.update(params)

plt.figure()

chirp_data[0,:,1] = chirp_data[0,:,2] + chirp_data[0,:,1]
chirp_data[0,:,0] = chirp_data[0,:,1] + chirp_data[0,:,0]
chirp_data[1,:,2] = chirp_data[0,:,0] + chirp_data[1,:,2]
chirp_data[1,:,1] = chirp_data[1,:,2] + chirp_data[1,:,1]
chirp_data[1,:,0] = chirp_data[1,:,1] + chirp_data[1,:,0]
chirp_data[2,:,2] = chirp_data[1,:,0] + chirp_data[2,:,2]
chirp_data[2,:,1] = chirp_data[2,:,2] + chirp_data[2,:,1]
chirp_data[2,:,0] = chirp_data[2,:,1] + chirp_data[2,:,0]
max_val = max(chirp_data[2,:,0])
chirp_data = chirp_data/max_val

ax = plt.gca()

ax.fill_between(chirp_vals, 0, chirp_data[2,:,0], color = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_data[2,:,1], color = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_data[2,:,1], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "////")
ax.fill_between(chirp_vals, 0, chirp_data[2,:,2], color = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_data[2,:,2], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "xxxx")

ax.fill_between(chirp_vals, 0, chirp_data[1,:,0], color = hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_data[1,:,1], color = hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_data[1,:,1], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "////")
ax.fill_between(chirp_vals, 0, chirp_data[1,:,2], color = hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_data[1,:,2], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "xxxx")

ax.fill_between(chirp_vals, 0, chirp_data[0,:,0], color = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_data[0,:,1], color = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_data[0,:,1], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "////")
ax.fill_between(chirp_vals, 0, chirp_data[0,:,2], color = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_data[0,:,2], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "xxxx")

plt.legend(proxies,labels,loc=1)

plt.gca().set_xlim(1.15,2.45)
plt.gca().set_xlabel("$\log\;M_{\mathrm{chirp}}\;\mathrm{[M_\odot]}$")
plt.gca().set_ylabel("relative number")

plt.savefig("merging_bh_formation.pdf", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                        transparent=False)

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
          'figure.subplot.right':0.95}

mpl.rcParams.update(params)

plt.figure()

chirp_LIGO_data[0,:,1] = chirp_LIGO_data[0,:,2] + chirp_LIGO_data[0,:,1]
chirp_LIGO_data[0,:,0] = chirp_LIGO_data[0,:,1] + chirp_LIGO_data[0,:,0]
chirp_LIGO_data[1,:,2] = chirp_LIGO_data[0,:,0] + chirp_LIGO_data[1,:,2]
chirp_LIGO_data[1,:,1] = chirp_LIGO_data[1,:,2] + chirp_LIGO_data[1,:,1]
chirp_LIGO_data[1,:,0] = chirp_LIGO_data[1,:,1] + chirp_LIGO_data[1,:,0]
chirp_LIGO_data[2,:,2] = chirp_LIGO_data[1,:,0] + chirp_LIGO_data[2,:,2]
chirp_LIGO_data[2,:,1] = chirp_LIGO_data[2,:,2] + chirp_LIGO_data[2,:,1]
chirp_LIGO_data[2,:,0] = chirp_LIGO_data[2,:,1] + chirp_LIGO_data[2,:,0]
max_val = max(chirp_LIGO_data[2,:,0])
chirp_LIGO_data = chirp_LIGO_data/max_val

ax = plt.gca()

ax.fill_between(chirp_vals, 0, chirp_LIGO_data[2,:,0], color = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data[2,:,1], color = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data[2,:,1], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "////")
ax.fill_between(chirp_vals, 0, chirp_LIGO_data[2,:,2], color = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data[2,:,2], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "xxxx")

ax.fill_between(chirp_vals, 0, chirp_LIGO_data[1,:,0], color = hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data[1,:,1], color = hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data[1,:,1], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "////")
ax.fill_between(chirp_vals, 0, chirp_LIGO_data[1,:,2], color = hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data[1,:,2], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "xxxx")

ax.fill_between(chirp_vals, 0, chirp_LIGO_data[0,:,0], color = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data[0,:,1], color = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data[0,:,1], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "////")
ax.fill_between(chirp_vals, 0, chirp_LIGO_data[0,:,2], color = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data[0,:,2], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "xxxx")

plt.legend(proxies,labels,loc=1)

plt.gca().set_xlim(1.15,2.45)
plt.gca().set_xlabel("$\log\;M_{\mathrm{chirp}}\;\mathrm{[M_\odot]}$")
plt.gca().set_ylabel("relative number")

plt.savefig("ligo_rates.pdf", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                        transparent=False)

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
          'figure.subplot.right':0.95}

mpl.rcParams.update(params)

plt.figure()

chirp_LIGO_data_low[0,:,1] = chirp_LIGO_data_low[0,:,2] + chirp_LIGO_data_low[0,:,1]
chirp_LIGO_data_low[0,:,0] = chirp_LIGO_data_low[0,:,1] + chirp_LIGO_data_low[0,:,0]
chirp_LIGO_data_low[1,:,2] = chirp_LIGO_data_low[0,:,0] + chirp_LIGO_data_low[1,:,2]
chirp_LIGO_data_low[1,:,1] = chirp_LIGO_data_low[1,:,2] + chirp_LIGO_data_low[1,:,1]
chirp_LIGO_data_low[1,:,0] = chirp_LIGO_data_low[1,:,1] + chirp_LIGO_data_low[1,:,0]
chirp_LIGO_data_low[2,:,2] = chirp_LIGO_data_low[1,:,0] + chirp_LIGO_data_low[2,:,2]
chirp_LIGO_data_low[2,:,1] = chirp_LIGO_data_low[2,:,2] + chirp_LIGO_data_low[2,:,1]
chirp_LIGO_data_low[2,:,0] = chirp_LIGO_data_low[2,:,1] + chirp_LIGO_data_low[2,:,0]
max_val = max(chirp_LIGO_data_low[2,:,0])
chirp_LIGO_data_low = chirp_LIGO_data_low/max_val

ax = plt.gca()

ax.fill_between(chirp_vals, 0, chirp_LIGO_data_low[2,:,0], color = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_low[2,:,1], color = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_low[2,:,1], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "////")
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_low[2,:,2], color = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_low[2,:,2], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "xxxx")

ax.fill_between(chirp_vals, 0, chirp_LIGO_data_low[1,:,0], color = hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_low[1,:,1], color = hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_low[1,:,1], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "////")
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_low[1,:,2], color = hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_low[1,:,2], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "xxxx")

ax.fill_between(chirp_vals, 0, chirp_LIGO_data_low[0,:,0], color = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_low[0,:,1], color = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_low[0,:,1], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "////")
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_low[0,:,2], color = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_low[0,:,2], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "xxxx")

plt.legend(proxies,labels,loc=1)

plt.gca().set_xlim(1.15,2.45)
plt.gca().set_xlabel("$\log\;M_{\mathrm{chirp}}\;\mathrm{[M_\odot]}$")
plt.gca().set_ylabel("relative number")

plt.savefig("ligo_rates_low.pdf", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                        transparent=False)

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
          'figure.subplot.right':0.95}

mpl.rcParams.update(params)

plt.figure()

chirp_LIGO_data_high[0,:,1] = chirp_LIGO_data_high[0,:,2] + chirp_LIGO_data_high[0,:,1]
chirp_LIGO_data_high[0,:,0] = chirp_LIGO_data_high[0,:,1] + chirp_LIGO_data_high[0,:,0]
chirp_LIGO_data_high[1,:,2] = chirp_LIGO_data_high[0,:,0] + chirp_LIGO_data_high[1,:,2]
chirp_LIGO_data_high[1,:,1] = chirp_LIGO_data_high[1,:,2] + chirp_LIGO_data_high[1,:,1]
chirp_LIGO_data_high[1,:,0] = chirp_LIGO_data_high[1,:,1] + chirp_LIGO_data_high[1,:,0]
chirp_LIGO_data_high[2,:,2] = chirp_LIGO_data_high[1,:,0] + chirp_LIGO_data_high[2,:,2]
chirp_LIGO_data_high[2,:,1] = chirp_LIGO_data_high[2,:,2] + chirp_LIGO_data_high[2,:,1]
chirp_LIGO_data_high[2,:,0] = chirp_LIGO_data_high[2,:,1] + chirp_LIGO_data_high[2,:,0]
max_val = max(chirp_LIGO_data_high[2,:,0])
chirp_LIGO_data_high = chirp_LIGO_data_high/max_val

ax = plt.gca()

ax.fill_between(chirp_vals, 0, chirp_LIGO_data_high[2,:,0], color = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_high[2,:,1], color = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_high[2,:,1], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "////")
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_high[2,:,2], color = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_high[2,:,2], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "xxxx")

ax.fill_between(chirp_vals, 0, chirp_LIGO_data_high[1,:,0], color = hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_high[1,:,1], color = hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_high[1,:,1], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "////")
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_high[1,:,2], color = hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_high[1,:,2], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "xxxx")

ax.fill_between(chirp_vals, 0, chirp_LIGO_data_high[0,:,0], color = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_high[0,:,1], color = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_high[0,:,1], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "////")
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_high[0,:,2], color = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2)
ax.fill_between(chirp_vals, 0, chirp_LIGO_data_high[0,:,2], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "xxxx")

plt.legend(proxies,labels,loc=1)

plt.gca().set_xlim(1.15,2.45)
plt.gca().set_xlabel("$\log\;M_{\mathrm{chirp}}\;\mathrm{[M_\odot]}$")
plt.gca().set_ylabel("relative number")

plt.savefig("ligo_rates_high.pdf", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                        transparent=False)

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
          'figure.subplot.right':0.95}

mpl.rcParams.update(params)

plt.figure()

merger_data[0,:,1] = merger_data[0,:,2] + merger_data[0,:,1]
merger_data[0,:,0] = merger_data[0,:,1] + merger_data[0,:,0]
merger_data[1,:,2] = merger_data[0,:,0] + merger_data[1,:,2]
merger_data[1,:,1] = merger_data[1,:,2] + merger_data[1,:,1]
merger_data[1,:,0] = merger_data[1,:,1] + merger_data[1,:,0]
merger_data[2,:,2] = merger_data[1,:,0] + merger_data[2,:,2]
merger_data[2,:,1] = merger_data[2,:,2] + merger_data[2,:,1]
merger_data[2,:,0] = merger_data[2,:,1] + merger_data[2,:,0]
max_val = max(merger_data[2,:,0])
merger_data = np.log10(merger_data/max_val+1e-99)

ax = plt.gca()

ax.fill_between(merger_vals, -100, merger_data[2,:,0], color = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2)
ax.fill_between(merger_vals, -100, merger_data[2,:,1], color = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2)
ax.fill_between(merger_vals, -100, merger_data[2,:,1], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "////")
ax.fill_between(merger_vals, -100, merger_data[2,:,2], color = hexcols[3], edgecolor =colorscale(hexcols[3],0.5), lw = 2)
ax.fill_between(merger_vals, -100, merger_data[2,:,2], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "xxxx")

ax.fill_between(merger_vals, -100, merger_data[1,:,0], color = hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2)
ax.fill_between(merger_vals, -100, merger_data[1,:,1], color = hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2)
ax.fill_between(merger_vals, -100, merger_data[1,:,1], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "////")
ax.fill_between(merger_vals, -100, merger_data[1,:,2], color = hexcols[5], edgecolor =colorscale(hexcols[5],0.5), lw = 2)
ax.fill_between(merger_vals, -100, merger_data[1,:,2], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "xxxx")

ax.fill_between(merger_vals, -100, merger_data[0,:,0], color = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2)
ax.fill_between(merger_vals, -100, merger_data[0,:,1], color = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2)
ax.fill_between(merger_vals, -100, merger_data[0,:,1], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "////")
ax.fill_between(merger_vals, -100, merger_data[0,:,2], color = hexcols[1], edgecolor =colorscale(hexcols[1],0.5), lw = 2)
ax.fill_between(merger_vals, -100, merger_data[0,:,2], color = "none", edgecolor = "k", linewidth = 0.0, hatch = "xxxx")

plt.legend(proxies,labels,loc=2)

plt.gca().set_xlim(-1,np.log10(13.8))
plt.gca().set_ylim(-2,0.3)

plt.gca().set_xlabel("$\log\;t_\mathrm{merger}\;\mathrm{[Gyr]}$")
plt.gca().set_ylabel("$\log$ (relative number)")

plt.savefig("tmerge_dist.pdf", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                        transparent=False)
