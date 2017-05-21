#!/usr/bin/env python
import os
import numpy as np
import sys
import shutil
import mesa_data as md
from mesa_data import *
from scipy import interpolate
import matplotlib.pyplot as plt
from pylab import *
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
          'figure.subplot.right':0.92}
mpl.rcParams.update(params)

hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',\
        '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466','#4477AA']

clight = 2.99792458e10
cgrav = 6.67428e-8
Msun = 1.9892e33
Rsun = 6.9598e10
Lsun = 3.8418e33
au = 1.495978921e13
secyear = 24*3600*365.25

def fRL(Md, Ma):
    q = Md/Ma
    return 0.49*np.power(q,2./3.)/(0.6*np.power(q,2./3.)+np.log(1+np.power(q,1./3.)))

#Peters (1964) beta, eq. (5.9)
def beta_gw(m1,m2):
    return 64./5.*cgrav**3/clight**5*m1*m2*(m1+m2)

#merger time (in Gyr) for initial orbital separation a (Rsun) and masses m1,m2 (Msun)
#Peters (1964), eq. (5.10)
def T_merger_a(a,m1,m2):
    return (a*Rsun)**4/(4.*beta_gw(m1,m2)*Msun**3)/(secyear*1e9)

lambda_names = ["lambda", "lambda_plus_Ecore", "lambda_sub_Uout"]
qratios = ["0.10","0.25","0.50","0.75","0.90"]
alpha_ces = ["0.20","0.60","1.00"]
files = {}
#These arrays store specific columns number to plot data at the end
colnums = []
colnums2 = []
colnums3 = []

for lambda_name in lambda_names:
    if not os.path.exists(lambda_name):
        os.makedirs(lambda_name)
    for alpha_ce_str in alpha_ces:
        if not os.path.exists(lambda_name+"/"+alpha_ce_str):
            os.makedirs(lambda_name+"/"+alpha_ce_str)
        files[lambda_name+alpha_ce_str] = open(lambda_name+"/"+alpha_ce_str+".dat","w")
        col_count = 0
        string =  '{0:>20}'.format("log_R(Rsun)"); col_count += 1
        colnums.append(col_count)
        string += '{0:>20}'.format("lambda_X=0.1"); col_count += 1
        string += '{0:>20}'.format("mass_X=0.1"); col_count += 1
        colnums.append(col_count)
        string += '{0:>20}'.format("lambda_X=0.2"); col_count += 1
        string += '{0:>20}'.format("mass_X=0.2"); col_count += 1
        colnums.append(col_count)
        string += '{0:>20}'.format("lambda_X=0.3"); col_count += 1
        string += '{0:>20}'.format("mass_X=0.3"); col_count += 1
        colnums.append(col_count)
        string += '{0:>20}'.format("lambda_relax"); col_count += 1
        string += '{0:>20}'.format("mass_relax"); col_count += 1
        string += '{0:>20}'.format("X_relax"); col_count += 1
        for qstr in qratios:
            colnums.append(col_count)
            string += '{0:>20}'.format("lambda_q="+qstr); col_count += 1
            string += '{0:>20}'.format("mass_q="+qstr); col_count += 1
            string += '{0:>20}'.format("X_q="+qstr); col_count += 1
            colnums3.append(col_count)
            string += '{0:>20}'.format("af_div_ai_q="+qstr); col_count += 1
            colnums2.append(col_count)
            string += '{0:>20}'.format("tmerge(Gyr)_q="+qstr); col_count += 1

        #print column numbers
        for i in range(col_count):
            files[lambda_name+alpha_ce_str].write('{0:>20}'.format(i+1))
        #print column names
        files[lambda_name+alpha_ce_str].write("\n"+string)

#initialize figures
fig, axes= plt.subplots(1)

models = sorted(os.listdir("models"))
for model in models:
    model_num = model[5:8]
    print("Checking model ",model_num)

    #check if simulation completed correctly
    outdata = os.popen('tail -n 100 simulations/'+model_num+'/out.txt').read()
    if outdata.find('Terminating because R>R_i/3') == -1:
        continue

    history = Mesa_Data("simulations/"+model_num+"/LOGS/history.data")

    surface_h1 = history.get("surface_h1")
    star_age = history.get("star_age")
    log_R = history.get("log_R")
    star_mass = history.get("star_mass")
    log_abs_mdot = history.get("log_abs_mdot")
    log_Lnuc = history.get("log_Lnuc")
    log_L = history.get("log_L")
    delta_Ecore = history.get("delta_Ecore")
    Ebind = history.get("Ebind")

    #indeces marking the switch from mass loss, to mass gain and thermal eq.
    k=0
    j=0
    #Find the switch to accretion
    for i, lmd in enumerate(history.get("log_abs_mdot")):
        if lmd < 0:
            k=i
            break
    #During accretion, get point at with thermal equilibrium is restored
    for i in range(k,len(history.get("log_abs_mdot"))):
        if abs(log_Lnuc[i]-log_L[i])<0.005 and star_age[i]-star_age[k]>1e5:
            j=i
            break

    #Find point where the thermally relaxed star would be just the size of the
    #one resulting from adiabatic mass loss
    masses_post = star_mass[j:]
    log_R_post = log_R[j:]
    f = interpolate.interp1d(masses_post, log_R_post)
    index_comp = -1
    for i in range(k,0,-1):
        if (star_mass[i]<min(masses_post) or star_mass[i]>max(masses_post)):
            continue
        if f(star_mass[i]) > log_R[i]:
            index_comp = i
            mass_comp = star_mass[i]
            X_comp = surface_h1[i]
            break

    for lambda_name in lambda_names:
        for alpha_ce_str in alpha_ces:
            fdata = files[lambda_name+alpha_ce_str]
            lambda_ce = history.get(lambda_name)
            string = '{0:>20e}'.format(float(model_num)/100)

            lambda_h3 = 0
            mass_h3 = 0
            lambda_h2 = 0
            mass_h2 = 0
            lambda_h1 = 0
            mass_h1 = 0
            #Get lambda values at fixed h1
            for i, h1 in enumerate(surface_h1):
                if mass_h3 == 0 and h1<0.3:
                    lambda_h3 = lambda_ce[i]
                    mass_h3 = star_mass[i]
                if mass_h2 == 0 and h1<0.2:
                    lambda_h2 = lambda_ce[i]
                    mass_h2 = star_mass[i]
                if mass_h1 == 0 and h1<0.1:
                    lambda_h1 = lambda_ce[i]
                    mass_h1 = star_mass[i]
                    break
            string += '{0:>20e}'.format(lambda_h1)
            string += '{0:>20e}'.format(mass_h1)
            string += '{0:>20e}'.format(lambda_h2)
            string += '{0:>20e}'.format(mass_h2)
            string += '{0:>20e}'.format(lambda_h3)
            string += '{0:>20e}'.format(mass_h3)

            #add values for the thermally relaxed boundary
            if index_comp >= 0:
                string += '{0:>20e}'.format(lambda_ce[index_comp])
                string += '{0:>20e}'.format(mass_comp)
                string += '{0:>20e}'.format(X_comp)
            else:
                string += '{0:>20e}'.format(0e0)
                string += '{0:>20e}'.format(0e0)
                string += '{0:>20e}'.format(0e0)

            axes.set_xlabel("$M\;[M_\\odot]$")
            axes.set_ylabel("$\\log~R/R_\\odot$")
            axes.plot(star_mass[0:k], log_R[0:k], label="Adiabatic")
            axes.plot(star_mass[j:], log_R[j:],":", label="Thermal eq.")

            alpha_ce = float(alpha_ce_str)
            #get lambda for inspiraling object
            masses_q = []
            for qstr in qratios:
                qratio = float(qstr)
                M2 = star_mass[0]*qratio
                Menv = star_mass[0]-star_mass
                rl = fRL(star_mass[0],M2)
                initial_separation = np.power(10,log_R[0])/rl
                final_separations = initial_separation*star_mass*M2/star_mass[0] \
                        /(M2+2*Menv/(alpha_ce*lambda_ce*rl))
                axes.plot(star_mass[0:k],np.log10(final_separations[0:k]*fRL(star_mass[0:k],M2)),"--",\
                        label = "q="+qstr, lw=1)
                min_separation = initial_separation
                for i, sep in enumerate(final_separations):
                    if sep>min_separation:
                        final_separations[i] = min_separation
                    else:
                        min_separation = sep
                frac = final_separations*fRL(star_mass,M2)/np.power(10,log_R)
                found_val = False
                for i in range(k):
                    if frac[i] > 1.01:
                        masses_q.append(star_mass[i])
                        string += '{0:>20e}'.format(lambda_ce[i])
                        string += '{0:>20e}'.format(star_mass[i])
                        string += '{0:>20e}'.format(surface_h1[i])
                        string += '{0:>20e}'.format(final_separations[i]/initial_separation)
                        string += '{0:>20e}'.format(T_merger_a(final_separations[i],star_mass[i],M2))
                        found_val = True
                        break
                if not found_val:
                    masses_q.append(0e0)
                    string += '{0:>20e}'.format(0e0)
                    string += '{0:>20e}'.format(0e0)
                    string += '{0:>20e}'.format(0e0)
                    string += '{0:>20e}'.format(0e0)
                    string += '{0:>20e}'.format(0e0)
            axes.legend(loc="lower right")
            axes.set_title("$\\log R_{\\rm i}/R_\\odot=$"+str(float(model_num)/100.)+" ($R="+str(int(np.power(10,float(model_num)/100.)))+"R_\\odot$)")
            plt.savefig(lambda_name+"/"+alpha_ce_str+"/R_"+model_num+".pdf")
            plt.cla()

            axes.set_xlabel("$M\;[M_\\odot]$")
            axes.set_ylabel("$\\log~E\;\\rm[erg]$")
            axes.plot([mass_h1,mass_h1],[-1000,1000],"--",color="0.8",label="X=0.1,0.2,0.3")
            axes.plot([mass_h2,mass_h2],[-1000,1000],"--",color="0.8")
            axes.plot([mass_h3,mass_h3],[-1000,1000],"--",color="0.8")
            axes.plot(star_mass[0:k], np.log10(-Ebind[0:k]), label="$-E_{\\rm bind}$")
            axes.plot(star_mass[0:k], np.log10(-Ebind[0:k]+delta_Ecore[0:k]),":", label="$-E_{\\rm bind}+\Delta E_{\\rm core}$")
            axes.plot([mass_comp,mass_comp],[-1000,1000],":",label="eq=ad")

            for i, qstr in enumerate(qratios):
                axes.plot([masses_q[i],masses_q[i]],[-1000,1000],"--",label="q="+qstr, lw=1)

            for i in range(k):
                if star_mass[i] < 0.99*star_mass[0]:
                    axes.set_ylim([max(np.log10(-Ebind[i]), np.log10(-Ebind[i]+delta_Ecore[i])),
                        max(np.log10(-Ebind[0:k])+0.1)])
                    break
            axes.set_xlim([star_mass[k],star_mass[0]])

            axes.set_title("$\\log R_{\\rm i}/R_\\odot=$"+str(float(model_num)/100.)+" ($R="+str(int(np.power(10,float(model_num)/100.)))+"R_\\odot$)")
            axes.legend(loc="upper right", fontsize=5, ncol=3)
            plt.savefig(lambda_name+"/"+alpha_ce_str+"/E_"+model_num+".pdf")
            plt.cla()

            #write one line of the table
            files[lambda_name+alpha_ce_str].write("\n"+string)


for lambda_name in lambda_names:
    for alpha_ce_str in alpha_ces:
        files[lambda_name+alpha_ce_str].close()
        data = np.genfromtxt(lambda_name+"/"+alpha_ce_str+".dat",skip_header=2, unpack=True)

        axes.set_xlabel("$R\;[R_\\odot]$")
        axes.set_ylabel("$\\lambda$")
        axes.plot(data[0], data[colnums[0]], label="$X=0.1$")
        axes.plot(data[0], data[colnums[1]], label="$X=0.2$")
        axes.plot(data[0], data[colnums[2]], label="$X=0.3$")
        axes.plot(data[0], data[colnums[3]], ":", label="$eq=ad$")
        num = 4
        for qstr in qratios:
            axes.plot(data[0], data[colnums[num]], "--", linewidth=1, label="$q="+qstr+"$")
            num += 1
        axes.set_title((lambda_name+" "+alpha_ce_str).replace("_","-"))
        axes.legend(loc="upper left", fontsize=5, ncol=3)
        axes.set_yscale("log", nonposy='clip')
        axes.set_ylim([0.005,50])
        plt.savefig(lambda_name+"/lambda_"+alpha_ce_str+".pdf")
        plt.cla()

        axes.set_xlabel("$R\;[R_\\odot]$")
        axes.set_ylabel("$t_{\\rm merge}\\;\\rm [Gyr]$")
        num = 0
        for qstr in qratios:
            axes.plot(data[0], data[colnums2[num]], "--", linewidth=1, label="$q="+qstr+"$")
            num += 1
        axes.set_title((lambda_name+" "+alpha_ce_str).replace("_","-"))
        axes.legend(loc="upper left", fontsize=5, ncol=3)
        axes.set_yscale("log", nonposy='clip')
        plt.savefig(lambda_name+"/tmerge_"+alpha_ce_str+".pdf")
        plt.cla()

        axes.set_xlabel("$R\;[R_\\odot]$")
        axes.set_ylabel("$a_{\\rm f}/a_{\\rm i}$")
        num = 0
        for qstr in qratios:
            axes.plot(data[0], data[colnums3[num]], "--", linewidth=1, label="$q="+qstr+"$")
            num += 1
        axes.set_title((lambda_name+" "+alpha_ce_str).replace("_","-"))
        axes.legend(loc="upper left", fontsize=5, ncol=3)
        axes.set_yscale("log", nonposy='clip')
        plt.savefig(lambda_name+"/af_div_ai_"+alpha_ce_str+".pdf")
        plt.cla()

