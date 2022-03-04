#We can compute Bulk flow and find the flow field inside at a radius r,given 
#the velocity field at the Surface(radius R),l is a 1-D array  

#l=0,1 mode not possible for Vlm_y and Vlm_psi,l=0 not possible for Vlm_phi.
def modify_(a,b,c,l): 
    x=np.where(l==1)
    y=np.where(l==0)
    for i in y:
        a[i]=complex(0,0)
        b[i]=complex(0,0)
        c[i]=complex(0,0)
    for i in x:
        a[i]=complex(0,0)
        b[i]=complex(0,0)
    
    return a,b,c

#To compute the bulk flow from the analytical results
def bulkflow(r_,Alm_,Clm_,l,R):
    #vslm, vtlm     = sht.analys(vtheta_s, vphi_s)
    #Computing the cofficients 
    #Alm            = vtlm/(R**l)
    #Clm            = (l*(l+1)* vslm)/(2*(R**(l+1)))
    r_lplus1       =r_**(l+1)                  #calculated r^(l+1)
    r_lminus1      =r_**(l-1)                  #Calculated r^(l-1)
    Vlm_phi        = Alm_*(r_**l)   
    Vlm_psi        = (Clm_*(-(R**2)*(r_lminus1)+(r_lplus1*(l+3)/(l+1)))/l)
    Vlm_y          = Clm_*(r_lplus1- (R**2)*r_lminus1)
    #the modes not allowed are removed , by making those cofficients zero.
    #Vlm_y_new,Vlm_psi_new,Vlm_phi_new=modify_(Vlm_y,Vlm_psi,Vlm_phi,l)
    vr,vtheta,vphi = sht.synth(Vlm_y,Vlm_psi,Vlm_phi)
    
    return  vr,vtheta,vphi


import numpy as np                                     
import shtns                                           
import cmath                                                
from mayavi.mlab import *
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import get_cmap, viridis
values = np.linspace(0., 1., 256)
lut = 255*get_cmap(viridis)(values.copy())
from mayavi import mlab

r_=np.linspace(0.2,1.,5)
lmax = 16
mmax=lmax
sht = shtns.sht(lmax, mmax)

nlats, nlons = sht.set_grid()

theta_vals = np.arccos(sht.cos_theta)
phi_vals   = (2.0*np.pi/nlons)*np.arange(nlons)
phi, theta = np.meshgrid(phi_vals, theta_vals)

l  = sht.l
m  = sht.m
l2 = l * (l + 1)

# given surface flow
vtheta_s      =np.sin(theta)#np.zeros_like(theta)#np.sin(theta)*np.sin(phi)#
vphi_s        =np.zeros_like(theta)#np.sin(theta)*np.cos(phi)##np.sin(theta)*np.cos(theta)
vslm, vtlm  = sht.analys(vtheta_s, vphi_s)
R=1
Alm = vtlm/(R**l)
Clm = (l*(l+1)* vslm)/(2*(R**(l+1)))
#phi, theta = np.meshgrid(
#  	np.linspace(0, 2 * np.pi, nlons ), 
#   	np.linspace(0, np.pi, nlats))

#Changing spherical coordinate field into Cartesian Velocity field

#now to plot the velocity field on the sphere
for i in range(len(r_)):
	ur,utheta,uphi=bulkflow(r_[i],Alm,Clm,l,R)
	x =r_[i]*np.sin(theta)*np.cos(phi)
	y =r_[i]*np.sin(theta)*np.sin(phi)
	z =r_[i]*np.cos(theta)
	ux =ur*np.sin(theta)*np.cos(phi)+utheta*np.cos(theta)*np.cos(phi)-uphi*np.sin(phi)
	uy =ur*np.sin(theta)*np.sin(phi)+utheta*np.cos(theta)*np.sin(phi)+uphi*np.cos(phi)                                             
	uz =ur*np.cos(theta)-utheta*np.sin(theta) 

	speed=np.sqrt(ux**2+uy**2+uz**2)
	
	#mlab.quiver3d(x, y, z,ux,uy,uz,colormap='jet')
	#mlab.vectorbar(title='Vector color mapping',orientation='vertical')
	
       
#mlab.show()
#def stress()

	 

	



