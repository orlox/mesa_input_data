#!/usr/bin/env python
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import matplotlib.patheffects as pe
from scipy.interpolate import spline
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

hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',\
        '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466','#4477AA']

Msun = 1.99e33
cgrav = 6.67e-8
Rsun = 6.96e10

#kepler third law to obtain P in days as function of a (in Rsun) and m1, m2 in Msun
def kepler3_P(a,m1,m2):
    return (4.*math.pi**2*(a*Rsun)**3/(cgrav*(m1+m2)*Msun))**(1./2.)/(24.*3600.)

def rlof_a(R,q):
    return R*(0.6*q**(-2./3.)+np.log(1+q**(-1./3.)))/(0.49*q**(-2./3.))

mpl.rcParams.update(params)

fig, axes= plt.subplots(1)

log_mass         = np.array([
        1.00,
        1.10,
        1.20,
#        1.30,
#        1.40,
#        1.50,
        1.55,
#        1.60,
        1.65,
        1.70,
        1.80,
        2.00,
        2.10
        ])
log_R  = np.array([
        4.2760535573793501E-001,
        4.9106730087380623E-001,
        5.5275024385163396E-001,
#        6.2462515401440166E-001,
#        6.8497589930410030E-001,
#        7.2680359096233449E-001,
        7.5452711142763518E-001,
#        8.0091387703154759E-001,
        8.0979203637750508E-001,
        8.3743230991846984E-001,
        8.9235613691518800E-001,
        1.0021433395949833E+000,
        1.0577850896107002E+000
        ])
omega_div_omegac = np.array([
        1.00,
        0.92,
        0.86,
#        0.78,
#        0.72,
#        0.66,
        0.58,
#        0.48,
        0.46,
        0.44,
        0.42,
        0.41,
        0.41
        ])
omegac = np.array([
        4.1382012380843602E-004,
        3.4093089160500161E-004,
        3.0480225177594429E-004,
#        2.8148196102212227E-004,
#        2.4658960024444325E-004,
#        2.1928813238945559E-004,
        2.0698775942728026E-004,
#        1.9617252597888118E-004,
        1.8445967748866469E-004,
        1.7402582772766377E-004,
        1.5394032876626964E-004,
        1.1800034838321288E-004,
        1.0232255503380916E-004
        ])

mass = np.linspace(10,100,100)
odoc = spline(log_mass,omega_div_omegac,np.log10(mass))
oc = spline(log_mass,omegac,np.log10(mass))
lr = spline(log_mass,log_R,np.log10(mass))

axes.fill_between(mass, 1, odoc, alpha=0.2)
#axes.fill_between(np.power(10,log_mass), 1, omega_div_omegac, alpha=0.5)
axes.plot(mass, 2*3.14/(0.7*24*3600)/oc,color=hexcols[6],label="Fixed $P_{\\rm rot,i}$")
axes.plot(mass, 2*3.14/(1.0*24*3600)/oc,color=hexcols[6])
axes.plot(mass, 2*3.14/(1.5*24*3600)/oc,color=hexcols[6])
axes.plot(mass,2*3.14/24/3600/oc/
        kepler3_P(rlof_a(np.power(10,lr),0.1),mass, mass*0.1),color=hexcols[2],ls="--",
        label="RLOF at ZAMS for fixed $q$"
        )
#axes.plot(mass,2*3.14/24/3600/oc/
#        kepler3_P(rlof_a(np.power(10,lr),0.25),mass, mass*0.25),color=hexcols[2],ls="--"
#        )
axes.plot(mass,2*3.14/24/3600/oc/
        kepler3_P(rlof_a(np.power(10,lr),0.9),mass, mass*0.9),color=hexcols[2],ls="--"
        )

axes.text(80, 0.47,'CHE',
     horizontalalignment='center',
     verticalalignment='top',
     fontsize=15,
     color="b", alpha=0.5, rotation=-2)

axes.text(80, 0.4,'Normal evolution',
     horizontalalignment='center',
     verticalalignment='top',
     fontsize=15,
     color="k", alpha=0.7, rotation=-2)

axes.text(30, 0.5,'$P_{\\rm rot,i}=0.7\;\\rm d$',
     horizontalalignment='center',
     verticalalignment='center',
     fontsize=10,
     color=hexcols[6], rotation=30)

axes.text(30, 0.36,'$P_{\\rm rot,i}=1.0\;\\rm d$',
     horizontalalignment='center',
     verticalalignment='center',
     fontsize=10,
     color=hexcols[6], rotation=20)

axes.text(30, 0.245,'$P_{\\rm rot,i}=1.5\;\\rm d$',
     horizontalalignment='center',
     verticalalignment='center',
     fontsize=10,
     color=hexcols[6], rotation=15)

axes.text(55, 0.695,'$q=0.1$',
     horizontalalignment='center',
     verticalalignment='center',
     fontsize=10,
     color=hexcols[2], rotation=10)

axes.text(55, 0.51,'$q=0.9$',
     horizontalalignment='center',
     verticalalignment='center',
     fontsize=10,
     color=hexcols[2], rotation=5)

axes.set_xlim([20,100])
axes.set_ylim([0.1,0.8])
axes.set_xlabel("$M_{\\rm i}\;[M_\odot]$")
axes.set_ylabel("$(\Omega/\Omega_{\\rm crit})_{\\rm surf,i}$")
axes.legend(loc="lower right", title = "$Z\\simeq Z_\odot/50$")

plt.savefig("../../images/che_diagram.pdf")

#        axarr[i].plot(profs_B[i].get("log_Teff")[3791:4273],profs_B[i].get("log_L")[3791:4273],color=hexcols[8],\
#                path_effects=[pe.Stroke(linewidth=7, foreground='k'), pe.Normal()], solid_capstyle='round',lw=6, zorder=-100)
