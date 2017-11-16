#!/usr/bin/env python3
import numpy as np
clight = 2.99792458e10
cgrav = 6.67428e-8
Msun = 1.9892e33
Rsun = 6.9598e10
Lsun = 3.8418e33
au = 1.495978921e13
secyear = 24*3600*365.25

x,y,z,am,hp,rho,vx,vy,vz,u,grpot,meanmolecular,cc, pressure = np.loadtxt('out0150.txt_20msun',unpack=True, skiprows=4)
h, he3, he4, c12, n14, o16, ne20, mg24 = np.loadtxt('compositions.sph_20msun',unpack=True)

x = x*Rsun
y = y*Rsun
z = z*Rsun
am = am*Msun
hp = hp*Rsun
rho = rho*Msun/Rsun**3
vx = vx*np.sqrt(cgrav*Msun/Rsun)
vy = vy*np.sqrt(cgrav*Msun/Rsun)
vz = vz*np.sqrt(cgrav*Msun/Rsun)
u = u*cgrav*Msun/Rsun
grpot = grpot*cgrav*Msun/Rsun
meanmolecular = meanmolecular/1.66e-24

#iterate to determine mass lost
total_mass_i = np.sum(am)
for i in range(100):
    total_E = 1./2.*am*(vx**2+vy**2+vz**2) + am*grpot + am*u
    #print(sum(am[total_E>0]), sum(am[total_E<0]/Msun))
    mask = total_E<0

    x_cm = np.sum(x[mask]*am[mask])/np.sum(am[mask])
    y_cm = np.sum(y[mask]*am[mask])/np.sum(am[mask])
    z_cm = np.sum(z[mask]*am[mask])/np.sum(am[mask])

    vx_cm = np.sum(vx[mask]*am[mask])/np.sum(am[mask])
    vy_cm = np.sum(vy[mask]*am[mask])/np.sum(am[mask])
    vz_cm = np.sum(vz[mask]*am[mask])/np.sum(am[mask])

    #displace everything to center of mass coords
    x = x-x_cm
    y = y-y_cm
    z = z-z_cm

    vx = vx-vx_cm
    vy = vy-vy_cm
    vz = vz-vz_cm

    total_mass = np.sum(am[mask])
    #print(total_mass_i/Msun,total_mass/Msun, vx_cm,vy_cm,vz_cm)

#only use particles that are bound
x = x[mask]
y = y[mask]
z = z[mask]
am = am[mask]
hp = hp[mask]
rho = rho[mask]
vx = vx[mask]
vy = vy[mask]
vz = vz[mask]
u = u[mask]
grpot = grpot[mask]
h = h[mask]
he3 = he3[mask]
he4 = he4[mask]
c12 = c12[mask]
n14 = n14[mask]
o16 = o16[mask]
ne20 = ne20[mask]
mg24 = mg24[mask]

#compute angular momentum
Jx = (y*vz-z*vy)*am
Jy = (z*vx-x*vz)*am
Jz = (x*vy-y*vx)*am

#project individual angular momenta into direction of total angular momentum
dirJ_x = np.sum(Jx)
dirJ_y = np.sum(Jy)
dirJ_z = np.sum(Jz)
norm = np.sqrt(dirJ_x**2+dirJ_y**2+dirJ_z**2)
dirJ_x = dirJ_x/norm
dirJ_y = dirJ_y/norm
dirJ_z = dirJ_z/norm
projected_J = Jx*dirJ_x + Jy*dirJ_y + Jz*dirJ_z

#compute moment of inertia for the z axis
Izz = (x*x+y*y)*am

#sort everything
r = np.sqrt(x**2+y**2+z**2)
order = np.argsort(r)
order = order[::-1]
x = x[order]
y = y[order]
z = z[order]
r = r[order]
am = am[order]
hp = hp[order]
rho = rho[order]
vx = vx[order]
vy = vy[order]
vz = vz[order]
u = u[order]
grpot = grpot[order]
h = h[order]
he3 = he3[order]
he4 = he4[order]
c12 = c12[order]
n14 = n14[order]
o16 = o16[order]
ne20 = ne20[order]
mg24 = mg24[order]
projected_J = projected_J[order]
Jx = Jx[order]
Jy = Jy[order]
Jz = Jz[order]
Izz = Izz[order]

vr = vx*x/r + vy*y/r + vz*z/r

#bin and output, include at least num_particles particles and dm mass per bin
num_particles = 100
dm = 0.2*Msun
bins = []
mass = 0
n = 0
for k in range(len(x)):
    if k == len(x)-1:
        bins[-1] = k
    else:
        mass = mass + am[k]
        n = n+1
        if n >= num_particles and mass >= dm:
            bins.append(k)
            mass = 0
            n = 0

outfile = open("entropy.dat","w")
outfile2 = open("composition.dat","w")
outfile3 = open("angular_momentum.dat","w")
#outfile.write(str(len(x)//num_particles+1)+"\n")
outfile.write(str(len(bins))+"\n")
outfile2.write(str(len(bins))+" 8\n")
outfile3.write(str(len(bins))+"\n")

rho_bin = np.zeros(len(bins))
r_bin = np.zeros(len(bins))
u_bin = np.zeros(len(bins))
#mass measured from the surface (M-M_r)
mass_bin = np.zeros(len(bins))
#these are specific angular momenta
j_bin = np.zeros(len(bins))
jx_bin = np.zeros(len(bins))
jy_bin = np.zeros(len(bins))
jz_bin = np.zeros(len(bins))
Izz_bin = np.zeros(len(bins))
vr_bin = np.zeros(len(bins))
#composition
h_bin = np.zeros(len(bins))
he3_bin = np.zeros(len(bins))
he4_bin = np.zeros(len(bins))
c12_bin = np.zeros(len(bins))
n14_bin = np.zeros(len(bins))
o16_bin = np.zeros(len(bins))
ne20_bin = np.zeros(len(bins))
mg24_bin = np.zeros(len(bins))
#Store arrays for all values except the central one.
#we obtain that one through extrapolation
mass = 0
for k in range(len(bins)-1):
    i = bins[k]
    j = bins[k+1]
    rho_bin[k] = np.sum(rho[i:j]*am[i:j])/np.sum(am[i:j])
    r_bin[k] = np.sum(r[i:j]*am[i:j])/np.sum(am[i:j])/Rsun
    u_bin[k] = np.sum(u[i:j]*am[i:j])/np.sum(am[i:j])
    j_bin[k] = np.sum(projected_J[i:j])/np.sum(am[i:j])
    jx_bin[k] = np.sum(Jx[i:j])/np.sum(am[i:j])
    jy_bin[k] = np.sum(Jy[i:j])/np.sum(am[i:j])
    jz_bin[k] = np.sum(Jz[i:j])/np.sum(am[i:j])
    Izz_bin[k] = np.sum(Izz[i:j])/np.sum(am[i:j])
    vr_bin[k] = np.sum(vr[i:j]*am[i:j])/np.sum(am[i:j])
    h_bin[k] = np.sum(h[i:j]*am[i:j])/np.sum(am[i:j])
    he3_bin[k] = np.sum(he3[i:j]*am[i:j])/np.sum(am[i:j])
    he4_bin[k] = np.sum(he4[i:j]*am[i:j])/np.sum(am[i:j])
    c12_bin[k] = np.sum(c12[i:j]*am[i:j])/np.sum(am[i:j])
    n14_bin[k] = np.sum(n14[i:j]*am[i:j])/np.sum(am[i:j])
    o16_bin[k] = np.sum(o16[i:j]*am[i:j])/np.sum(am[i:j])
    ne20_bin[k] = np.sum(ne20[i:j]*am[i:j])/np.sum(am[i:j])
    mg24_bin[k] = np.sum(mg24[i:j]*am[i:j])/np.sum(am[i:j])
    mass_bin[k] = mass
    mass = mass + np.sum(am[i:j])/Msun

    i = i + num_particles

mass_bin[-1] = total_mass/Msun

#get sizes of bins
dq = np.zeros(len(bins))
for k in range(len(bins)-1):
    dq[k] = abs(mass_bin[k]-mass_bin[k+1])/mass_bin[-1]

#get face centered values
dx = np.zeros(len(bins))
dx[0] = 0.5*dq[0]
print(dq[0],dx[0])
for k in range(1,len(bins)-1):
    dx[k] = dx[k-1] + 0.5*(dq[k] + dq[k-1])
    #print(dq[k],dx[k])
dx[-1] = 1

#extrapolate internal energy
u_bin[-1] = u_bin[-2] + (dx[-1]-dx[-2])*(u_bin[-2]-u_bin[-4])/(dx[-2]-dx[-4])
#extrapolate density (in log space)
log_rho_bin = np.log10(rho_bin)
log_rho_bin[-1] = log_rho_bin[-2] + (dx[-1]-dx[-2])*(log_rho_bin[-2]-log_rho_bin[-4])/(dx[-2]-dx[-4])
rho_bin[-1] = np.power(10,log_rho_bin[-1])
#"extrapolate" composition
h_bin[-1] = h_bin[-2]
he3_bin[-1] = he3_bin[-2]
he4_bin[-1] = he4_bin[-2]
c12_bin[-1] = c12_bin[-2]
n14_bin[-1] = n14_bin[-2]
o16_bin[-1] = o16_bin[-2]
ne20_bin[-1] = ne20_bin[-2]
mg24_bin[-1] = mg24_bin[-2]
#central angular momentum must be zero
j_bin[-1] = 0
jx_bin[-1] = 0
jy_bin[-1] = 0
jz_bin[-1] = 0
Izz_bin[-1] = 0
vr_bin[-1] = 0


for k in range(len(bins)):
    outfile.write(str(dx[k])+" "+str(rho_bin[k]) + " " + str(u_bin[k]))
    outfile2.write('{0:20.10}'.format(dx[k])+" "+'{0:20.10}'.format(h_bin[k])+" "+'{0:20.10}'.format(he3_bin[k])+
            " "+'{0:20.10}'.format(he4_bin[k])+" "+'{0:20.10}'.format(c12_bin[k])+
            " "+'{0:20.10}'.format(n14_bin[k])+" "+'{0:20.10}'.format(o16_bin[k])+
            " "+'{0:20.10}'.format(ne20_bin[k])+" "+'{0:20.10}'.format(mg24_bin[k]))
    outfile3.write(str(dx[k])+ " " + str(j_bin[k]))
    outfile.write("\n")
    outfile2.write("\n")
    outfile3.write("\n")

    print(k,total_mass/Msun-mass_bin[k],r_bin[k],rho_bin[k],u_bin[k],\
            j_bin[k],jx_bin[k],jy_bin[k],jz_bin[k], Izz_bin[k], j_bin[k]/Izz_bin[k], np.sqrt(cgrav*(total_mass-Msun*mass_bin[k])/(r_bin[k]*Rsun)**3),\
            h_bin[k],he4_bin[k],c12_bin[k]+n14_bin[k]+o16_bin[k]+ne20_bin[k]+mg24_bin[k], vr_bin[k])

outfile.close()

