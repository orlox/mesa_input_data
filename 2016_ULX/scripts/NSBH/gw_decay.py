import numpy as np
from scipy import integrate
from constants import *
import math

#kepler third law to obtain a in Rsun as function of P (in days) and m1, m2 in Msun
def kepler3_a(P,m1,m2):
    return float(((P*24.*3600.)**2*cgrav*(m1+m2)*Msun/(4.*math.pi**2))**(1./3.)/Rsun);

#kepler third law to obtain P in days as function of a (in Rsun) and m1, m2 in Msun
def kepler3_P(a,m1,m2):
    return float((4.*math.pi**2*(a*Rsun)**3/(cgrav*(m1+m2)*Msun))**(1./2.)/(24.*3600.));

#Peters (1964) beta, eq. (5.9)
def beta_gw(m1,m2):
    return 64./5.*cgrav**3/clight**5*m1*m2*(m1+m2)

#Peters (1964) c0_gw, eq. (5.11)
def c0_gw(a0,e0):
    return a0*(1-e0**2)/(e0**(12./19.)*(1+121./304.*e0**2)**(870./2299.))

#merger time (in Gyr) for initial orbital separation a (Rsun) and masses m1,m2 (Msun)
#Peters (1964), eq. (5.10)
def T_merger_a(a,m1,m2):
    return (a*Rsun)**4/(4.*beta_gw(m1,m2)*Msun**3)/(secyear*1e9)

#merger time (in Gyr) for initial orbital period P (days) and masses m1,m2 (Msun)
def T_merger_P(P,m1,m2):
    return T_merger_a(kepler3_a(P,m1,m2),m1,m2)

#integrand in eq. (5.14) of Peters (1964)
T_merger_integrand = lambda e: e**(29./19.)*(1+121./304.*e**2)**(1181./2299.)/(1-e**2)**1.5

#Merger time (in Gyr) for initial orbital separation a (Rsun) and eccentricity e.
#This integrates de/dt and da/dt
#Timestep is defined as dt_control*min(a/(da/dt),e/(de/dt)).
#Evolution is terminated once the separation is equal to end_condition*max(R_1sch,R_2sch)
#Timestep limit for eccentricity is only used while e>emin
def T_merger_a_e(a,e,m1,m2):
    #if input eccentricity is zero just return the exact result
    if e == 0:
        return T_merger_a(a,m1,m2)

    a = a*Rsun
    m1 = m1*Msun
    m2 = m2*Msun
    
    #integral can be very small for small eccentricity, so fix tolerance
    #based on the ratio of c0_gw(a,e)**4/(a)**4
    c0 = c0_gw(a,e)
    T = 12./19.*c0**4/beta_gw(m1,m2)*\
            integrate.quadrature(T_merger_integrand,0,e,tol=1e-10/(c0**4/a**4))[0]

    return T/(secyear*1e9)

#Merger time (in Gyr) for initial orbital period P (days) and eccentricity e.
def T_merger_P_e(P,e,m1,m2):
    return T_merger_a_e(kepler3_a(P,m1,m2),e,m1,m2)

