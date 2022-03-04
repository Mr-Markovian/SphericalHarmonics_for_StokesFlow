# Bulk flow at a particular r,l is a 1-D array  
def modify_(a,b,c,l): #l=0,1 mode not possible for Vlm_y and Vlm_psi,l=0 not possible for Vlm_phi.
    x=np.where(l==0)
    for i in x:
        np.append(complex(0,0),a)
        np.append(complex(0,0),b)
 	np.append(complex(0,0),c)   
    return a,b,c

def bulkflow(r_,vtheta_s,vphi_s,l,R):
    vslm, vtlm     = sht.analys(vtheta_s, vphi_s)
    lp=l[1:]
    #Computing the cofficients 
    #we can't divide by 0
    Alm            = vtlm[1:]/(R**lp)
    Clm            = (l2[1:]* vslm[1:])/(2*(R**(lp+1)))
    r_lplus1       =r_**(lp+1)                  #calculated r^(l+1)
    r_lminus1      =r_**(lp-1)                  #Calculated r^(l-1)
    Vlm_phi        = Alm*(r_**lp)   
    Vlm_psi        = (Clm*(-(R**2)*(r_lminus1)+(r_lplus1*(lp+3)/(lp+1)))/lp)
    Vlm_y          = Clm*(r_lplus1- (R**2)*r_lminus1)
    #the l=0 modes not allowed are removed , by making those cofficients zero.
    Vlm_y_new,Vlm_psi_new,Vlm_phi_new=modify_(Vlm_y,Vlm_psi,Vlm_phi,l)
    #vr,vtheta,vphi =sht.synth(Vlm_y,Vlm_psi,Vlm_phi)
    
    return  Vlm_y_new,Vlm_psi_new,Vlm_phi_new  


import numpy as np                                     
import shtns                                           
import cmath                                                
from mayavi.mlab import *
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import get_cmap, viridis
values = np.linspace(0., 1., 256)
lut = 255*get_cmap(viridis)(values.copy())
from mayavi import mlab
import matplotlib.pyplot as plt
r=np.linspace(0.2,1.,5)
lmax = 50
mmax=lmax
sht = shtns.sht(lmax, mmax)

nlats, nlons = sht.set_grid()

theta_vals = np.arccos(sht.cos_theta)
phi_vals   = (2.0*np.pi/nlons)*(np.arange(nlons)-nlons/2.)
phi, theta = np.meshgrid(phi_vals, theta_vals)

l = sht.l
m = sht.m
l2 = l * (l + 1)

# Given surface flow
vtheta = np.zeros_like(phi)
vphi = 3.0 * np.cos(theta)

Vslm, Vtlm = sht.analys(vtheta, vphi)
utheta,uphi= sht.synth(Vslm,Vtlm) #ur, utheta , uphi
#ur,utheta,uphi=sht.synth(urlm,upsi_lm,uphi_lm)
#vr_,vtheta_,vphi_=sht.synth(urlm,Vslm,Vtlm)
vslm,vtlm=sht.analys(utheta,uphi)
vtheta_,vphi_= sht.synth(vslm,vtlm)
print(vtheta-vtheta_)
print(vphi-vphi_)
#calculating the spectral components of the difference itself
dvslm,dvtlm=sht.analys(vtheta-utheta,vphi-uphi)
#print(dvslm)
#print(dvtlm)
dvtheta,dvphi=sht.synth(dvslm,dvtlm)
#print(dvtheta)
#Changing spherical coordinate field into Cartesian Velocity field
#vx=ur*np.sin(theta)*np.cos(phi)+vtheta*np.cos(theta)*np.cos(phi)-vphi*np.sin(phi)
#vy=ur*np.sin(theta)*np.sin(phi)+vtheta*np.cos(theta)*np.sin(phi)+vphi*np.cos(phi)                                                  
#vz=ur*np.cos(theta)-utheta*np.sin(theta)  

#now to plot the velocity field on the sphere

#x=np.sin(theta)*np.cos(phi)
#y=np.sin(theta)*np.sin(phi)
#z=np.cos(theta)

#mlab.quiver3d(x, y, z,vx,vy,vz)
#mlab.show()




