from gw_decay import *
from constants import *
import numpy as np
import math
from scipy.stats import maxwell

#Final orbital parameters in a circular sistem with a separation a (Rsun), masses m1, m2 (Msun),
#where m1 ejects mass and ends up with a mass m1f, and a kick velocity w (km/s). This ignores
#impulse velocities. Kick is in direction theta, phi in radians.
#Taken from Tauris+ (1999), MNRAS, 310, 1165-1169
def post_kick_parameters_a(a,m1,m2,m1f,w,theta,phi):
    m1 = m1*Msun
    m2 = m2*Msun
    m1f = m1f*Msun
    a = a*Rsun
    w = w*1e5
    mtilde = (m1f+m2)/(m1+m2)
    v = np.sqrt(cgrav*(m1+m2)/a)

    xi = (v**2+w**2+2*w*v*np.cos(theta))/(mtilde*v**2)
    Q = xi - 1 - (w*np.sin(theta)*np.cos(phi))**2/(mtilde*v**2)

    afinal = a/(2-xi)
    efinal = np.sqrt(1+(xi-2)*(Q+1))

    return [afinal/Rsun,efinal]

#Merger time (in Gyr) for initial orbital period P (days) and eccentricity e.
def post_kick_parameters_P(P,m1,m2,m1f,w,theta,phi):
    orbit = post_kick_parameters_a(kepler3_a(P,m1,m2),m1,m2,m1f,w,theta,phi)
    return [kepler3_P(orbit[0],m1f,m2), orbit[1]]

def flat_vdist(x):
    return x*300.

def hobbs_vdist(x):
    return maxwell.ppf(x)*265.

def hobbs_vdist2(x):
    return maxwell.ppf(x)*265./10.

#sample a velocity distribution vdist and return merger times, final periods and eccentricities
#vdist is the inverse cumulative distribution velocity
def sample_kick_distribution_a(a,m1,m2,m1f,vdist=hobbs_vdist,\
        num_v=100,num_theta=100,num_phi=20,filename="kick.data",\
        report_progress=False, save_data=False, output_in_days=False):
    linspace_v = np.linspace(0,1,num_v+2)
    v_vals = vdist(linspace_v[1:-1])

    theta_vals = np.linspace(0,math.pi,num_theta+2)[1:-1]
    phi_vals = np.linspace(0,2*math.pi,num_phi+2)[1:-1]

    dtheta = math.pi/(num_theta)
    dphi = 2*math.pi/(num_phi)

    if save_data:
        f = open(filename, 'w')

    k = 0
    numsim = len(v_vals)*len(theta_vals)*len(phi_vals)
    sum_weights = 0
    sum_weights_merge = 0
    sum_weights_disrupt = 0
    for i, v_val in enumerate(v_vals):
        for theta_val in theta_vals:
            for phi_val in phi_vals:
                k=k+1
                if report_progress and k%10000==0:
                    print(str(k)+"/"+str(numsim)+", "+str(v_val)+", "+str(theta_val)+", "+str(phi_val))
                weight = dtheta*dphi*np.sin(theta_val)/(4*math.pi)/(num_v)
                sum_weights += weight
                orbit = post_kick_parameters_a(a,m1,m2,m1f,v_val,theta_val,phi_val)
                if orbit[1] > 1 or not np.isfinite(orbit[0]) or not np.isfinite(orbit[1]):
                    sum_weights_disrupt += weight
                    merge_time = float("inf")
                else:
                    merge_time = T_merger_a_e(orbit[0],orbit[1],m1f,m2)

                if output_in_days:
                    orbit[0] = kepler3_P(orbit[0],m1f,m2)
                if save_data:
                    f.write("{0:>20e}".format(v_val)+"   "\
                           +"{0:>20e}".format(theta_val)+"   "\
                           +"{0:>20e}".format(phi_val)+"   "\
                           +"{0:>20e}".format(orbit[0])+"   "\
                           +"{0:>20e}".format(orbit[1])+"   "\
                           +"{0:>20e}".format(merge_time)+"   "\
                           +"{0:>20e}".format(weight))
                    f.write("\n")

                if merge_time < 13.8:
                    sum_weights_merge += weight

    if save_data:
        f.close()

    if report_progress:
        if not output_in_days:
            print("modelled system with a=",a,"m1=",m1,"m2=",m2,"m1f=",m1f)
        else:
            print("modelled system with P=",kepler3_P(a,m1,m2),"m1=",m1,"m2=",m2,"m1f=",m1f)
        print("RESULTS:")
        print("sum weights:",sum_weights)
        print("merging fraction:",sum_weights_merge/sum_weights)
        print("disrupt fraction:",sum_weights_disrupt/sum_weights)
    return [sum_weights_merge/sum_weights, sum_weights_disrupt/sum_weights]

def sample_kick_distribution_P(P,m1,m2,m1f,vdist=hobbs_vdist,\
        num_v=200,num_theta=200,num_phi=20,filename="kick.data",\
        report_progress=False, save_data=False):
    return sample_kick_distribution_a(kepler3_a(P,m1,m2),m1,m2,m1f,vdist=vdist,\
            num_v=num_v,num_theta=num_theta,num_phi=num_phi,filename=filename,\
            report_progress=report_progress, save_data=save_data, output_in_days=True)
