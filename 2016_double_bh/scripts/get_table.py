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

def format_entry(data):
    str = ""
    for datum in data:
        if type(datum) == type("") or type(datum) == type(1):
            str += '{0:>20}'.format(datum)
        else:
            str += '{0:>20e}'.format(datum)
    return str

def generate_data(mass_donor, period, qratio, metallicity, result,\
        had_contact = "-", M_1f = "-", M_2f = "-", P_f = "-", KP_1 = "-", KP_2 = "-", merge_time = "-"):
    return {
            "log(M_1i)":mass_donor,
            "qratio":qratio,
            "P_i":period,
            "Z":float(metallicity),
            "result":result,
            "had_contact":had_contact,
            "M_1f":M_1f,
            "M_2f":M_2f,
            "P_f":P_f,
            "Kerr_param_1":KP_1,
            "Kerr_param_2":KP_2,
            "merge_time":merge_time
            }

def create_table(masses_donor, periods, qratio, metallicity, folder ='.', table_name = "data.dat"):
    folder = folder + "/"
    data_file = open(table_name, 'w')

    data_file.write("#Data from case M models\n")
    data_file.write("#Possible outcomes are contained in the 'result' column and are:\n")
    data_file.write("#    -ZAMS_L2_overflow: Contact at ZAMS would be deep enough to reach L2. Should merge.\n")
    data_file.write("#    -L2_overflow: Stars undergo deep contact reaching L2. Should merge.\n")
    data_file.write("#    -off_CHE: Rotational mixing is not efficient enough and one component evolves not chemichally homogeneous.\n")
    data_file.write("#    -PISN: One star in the system depletes helium inside the PISN range (60-130 Msun).\n")
    data_file.write("#    -double_BH: System should form double BHs (both stars deplete helium outside 60-130 Msun PISN range).\n")
    data_file.write("#    -convergence_error: Run did not complete due to convergence issues.\n")
    data_file.write("#    -missing_data: For misterious reasons the data for this model is missing/incomplete.\n")
    data_file.write("#The column had_contact specifies whether the system had a contact phase, it can have 3 values:\n")
    data_file.write("#    -0: No contact.\n")
    data_file.write("#    -1: Reached contact during the main sequence.\n")
    data_file.write("#    -2: System already at contact during ZAMS.\n")
    data_file.write("#Merge times and kerr parameters are provided assuming direct collapse. Even when a PISN is expected,\n")
    data_file.write("#these numbers are given due to possible uncertainties in the range for occurrence of PISN.\n")
    data_file.write(format_entry(["log10(M_1i)(Msun)","qratio(M_2i/M_1i)","P_i(days)","metallicity","result","had_contact","M_1f(Msun)","M_2f(Msun)","P_f(days)","merge_time(Gyr)","Kerr_param_1","Kerr_param_2"]))
    data_file.write("\n")

    for i, period in enumerate(periods):
        for j, mass_donor in enumerate(masses_donor):
            result = ""
            data = None
            had_contact = 0
            folder_name = folder+mass_donor+"_"+qratio+"_"+period+"/"
            if not os.path.isdir(folder_name):
                result = "missing_data"
                data = generate_data(mass_donor, period, qratio, metallicity, result)
            else:
                outdata = os.popen('tail -n 40 '+folder_name+'out.txt').read()
                outfile = open(folder_name+'out.txt')
                s = mmap.mmap(outfile.fileno(), 0, access=mmap.ACCESS_READ)
                if outdata.find('DATE:') == -1:
                    result = "missing_data"
                    data = generate_data(mass_donor, period, qratio, metallicity, result)
                else:
                    if s.rfind('model is overflowing at ZAMS') != -1:
                        had_contact = 2
                    if outdata.find('Terminate due to overflow of L2 at ZAMS') != -1 or \
                            float(period)<0.9 and float(mass_donor)>2.3:
                        result = "ZAMS_L2_overflow"
                        data = generate_data(mass_donor, period, qratio, metallicity, result, had_contact = 2)
                    elif outdata.find('Terminate due to secondary not evolving homogeneously') != -1:
                        result = "off_CHE"
                        data = generate_data(mass_donor, period, qratio, metallicity, result)
                    elif outdata.find('Terminate due to primary not evolving homogeneously') != -1:
                        result = "off_CHE"
                        data = generate_data(mass_donor, period, qratio, metallicity, result)
                    elif outdata.find('termination code: Terminate because of L2 overflow') != -1:
                        result = "L2_overflow"
                        if had_contact == 0:
                            had_contact = 2
                        data = generate_data(mass_donor, period, qratio, metallicity, result, had_contact = had_contact)
                    elif outdata.find('Terminate due to helium depletion') != -1:
                        binary_history = history_data(folder_name, slname = "binary_history.data", clean_starlog=False)
                        star_history1 = history_data(folder_name+"/LOGS1", slname = "history.data", clean_starlog=False)
                        star_history2 = history_data(folder_name+"/LOGS2", slname = "history.data", clean_starlog=False)
                        mass_end1 = binary_history.get("star_1_mass")[-1]
                        mass_end2 = binary_history.get("star_2_mass")[-1]
                        period_end = binary_history.get("period_days")[-1]
                        kerr_param1 = np.power(10,star_history1.get("log_total_angular_momentum")[-1])\
                                *clight/((mass_end1*msun)**2*cgrav)
                        kerr_param2 = np.power(10,star_history2.get("log_total_angular_momentum")[-2])\
                                *clight/((mass_end2*msun)**2*cgrav)

                        #merge_time2 = pow(10, 8.0/3.0*log10(period_end)+2.276-5.0/3.0*log10(mass_end1+mass_end2))
                        merge_time = 189 * pow(period_end,8.0/3.0)*pow(mass_end1+mass_end2,-5.0/3.0)

                        if had_contact == 0:
                            overflow_1 = binary_history.get("rl_relative_overflow_1")
                            overflow_2 = binary_history.get("rl_relative_overflow_2")
                            lg_mstar_dot_1 = binary_history.get("lg_mstar_dot_1")
                            for k in range(len(overflow_1)):
                                if overflow_1[k]>0 and overflow_2[k]>0 and lg_mstar_dot_1[k]>-90:
                                    had_contact = 1
                                    break

                        if (mass_end1 > 60 and mass_end1 < 130) or (mass_end2 > 60 and mass_end2 < 130):
                            result = "PISN"
                        else:
                            result = "double_BH"

                        data =  generate_data(mass_donor, period, qratio, metallicity, result,\
                            had_contact = had_contact, M_1f = mass_end1, M_2f = mass_end2, P_f = period_end,\
                            KP_1 = kerr_param1, KP_2 = kerr_param2, merge_time = merge_time)

                    else:
                        result = "convergence_error"
                        data = generate_data(mass_donor, period, qratio, metallicity, result)

                print folder_name
                data_file.write(format_entry([data["log(M_1i)"],data["qratio"],data["P_i"],data["Z"],\
                        data["result"],data["had_contact"],data["M_1f"],data["M_2f"],data["P_f"],\
                        data["merge_time"],data["Kerr_param_1"],data["Kerr_param_2"]]))
                data_file.write("\n")
    data_file.close()



lg_masses_donor =    ["1.4","1.5","1.6","1.7","1.8","1.9","2.0","2.1","2.2","2.3","2.4","2.5","2.6","2.7"]
periods = ["0.50","0.60","0.70","0.80","0.90","1.00","1.10","1.20"\
        ,"1.30","1.40","1.50","1.60","1.70","1.80","1.90"\
        ,"2.00","2.10","2.20","2.30","2.40","2.50","2.60","2.70","2.80","2.90","3.00"]

qratios = ["1.0","0.9","0.8"]

folders = ["Z50","Z20","Z10"]
folders = ["Z50"]
metallicities = {
        "Z50" : 0.017/50,
        "Z20" : 0.017/20,
        "Z10" : 0.017/10
        }

for folder in folders:
    for qratio in qratios:
        folder_name = "double_bh_"+ folder
        create_table(lg_masses_donor, periods, qratio, metallicities[folder], folder = folder_name, table_name = "data_"+folder+"_"+qratio+".dat")
