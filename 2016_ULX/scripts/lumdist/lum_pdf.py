#!/usr/bin/env python
import numpy as np
import math
import scipy
import itertools
from scipy.stats import poisson

clight = 3e10
cgrav = 6.67e-8
msun = 2e33
Lsun = 3.9e33


folders = ["-2.500","-3.000","-3.500","-4.000"]#,"-4.500","-5.000","-5.500","-6.000"]
folders = ["-2.500","-3.000","-4.500","-5.000"]
qratios = ["0.050","0.100","0.150","0.200","0.250",\
        "0.300","0.350","0.400","0.450","0.500","0.550","0.600"]
for folder in folders:
    datafile = open("plots_lum/lumdist_"+folder+".dat","w")
    datafile.write("{0:>20}".format("SFR")+"   "+\
            "{0:>20}".format("median_Lx")+"   "+\
            "{0:>20}".format("expected_Lx")+"   "+\
            "{0:>20}".format("perc25_Lx")+"   "+\
            "{0:>20}".format("perc75_Lx")+"   "+\
            "{0:>20}".format("perc10_Lx")+"   "+\
            "{0:>20}".format("perc20_Lx")+"   "+\
            "{0:>20}".format("perc30_Lx")+"   "+\
            "{0:>20}".format("perc40_Lx")+"   "+\
            "{0:>20}".format("perc50_Lx")+"   "+\
            "{0:>20}".format("perc60_Lx")+"   "+\
            "{0:>20}".format("perc70_Lx")+"   "+\
            "{0:>20}".format("perc80_Lx")+"   "+\
            "{0:>20}".format("perc90_Lx")
            )
    datafile.write("\n")
    luminosities = []
    bh_masses = []
    weights = []
    for qratio in qratios:
        file_name = "mt_data_"+folder+"_"+qratio+".dat"
        print file_name
        try:
            data = np.loadtxt(file_name, skiprows=2, unpack = True)
            if data.size == 0: continue
        except:
            print "cant open file, or no data", file_name
            continue
        luminosities.append(data[11])
        bh_masses.append(data[6])
        weights.append(data[4]*data[8])

    luminosities = np.power(10,np.concatenate(luminosities))*Lsun/1e39
    bh_masses = np.log10(np.concatenate(bh_masses))
    weights = np.concatenate(weights)

    SFRs = np.power(10,np.linspace(-1,3,161))
    mean_luminosity_per_sfr = sum(luminosities*weights*0.01/3.0)
    for k, SFR in enumerate(SFRs):
        Nsim = 10000
        int_L = np.zeros(Nsim)
        print "modeling SFR", SFR, np.log10(SFR)
        for i in range(len(luminosities)):
            int_L += poisson.rvs(weights[i]*SFR*0.01/3, size=Nsim)*luminosities[i]
        print SFR, np.median(int_L), np.percentile(int_L,10), np.percentile(int_L,90), SFR*mean_luminosity_per_sfr
        datafile.write("{0:>20e}".format(SFR)+"   "+\
                "{0:>20e}".format(SFR*mean_luminosity_per_sfr)+"   "+\
                "{0:>20e}".format(np.median(int_L))+"   "+\
                "{0:>20e}".format(np.percentile(int_L,25))+"   "+\
                "{0:>20e}".format(np.percentile(int_L,75))+"   "+\
                "{0:>20e}".format(np.percentile(int_L,10))+"   "+\
                "{0:>20e}".format(np.percentile(int_L,20))+"   "+\
                "{0:>20e}".format(np.percentile(int_L,30))+"   "+\
                "{0:>20e}".format(np.percentile(int_L,40))+"   "+\
                "{0:>20e}".format(np.percentile(int_L,50))+"   "+\
                "{0:>20e}".format(np.percentile(int_L,60))+"   "+\
                "{0:>20e}".format(np.percentile(int_L,70))+"   "+\
                "{0:>20e}".format(np.percentile(int_L,80))+"   "+\
                "{0:>20e}".format(np.percentile(int_L,90))
                )
        datafile.write("\n")
        datafile.flush()

    datafile.close()
