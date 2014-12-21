# -*- coding: utf-8 -*-
"""
Created on Sat Dec 20 13:53:32 2014

@author: michael rottmaier
"""

import numpy                       
import matplotlib.pyplot as plt   

def V(rho, V_max, rho_max):
    """Calculates the velocity;
    
    Returns the velocity V = V_max*(1 - rho/rho_max)  where 
    * V_max is the max velocity  
    * rho is the traffic density 
    * rho_max is the maximal traffic density
    
    
    Parameters
    ----------
    
    rho : float
        traffic density
    rho_max : float
        max traffic density
    V_max   : float 
        max velocity
    Returns
    -------
    
    V : float 
    """     
    return V_max*(1 - rho/rho_max)


def traffic(T, V_max, L, rho_max, rho0, boundary_0, nx):
    """Solve the traffic flow equation.
    
    Solves the equation d_t rho + (d_rho F) d_x rho = 0 where 
    * F(rho) is the traffic Flux
    * V_max, L, rho_max are model parameters 
    * the domain is x \in [0, L]
    * $T/(60*\Delta t) timesteps are taken, with \Delta t = 0.001 [hour]
    * the initial data is the hat function
    
    Parameters
    ----------
    
    T : integer
        Time [min] to simulate 
        
    Returns
    -------
    
    rho : traffic density  
    """ 
    
    #Initialize           
    dx = L/(nx-1)
    dt = .001                   #timestep [hour]
    nt = int(T/(60*dt))+1
   

    
    #-----finite difference method ----------
    rho = rho0.copy()
    
    for n in range(nt): 
        rhon = rho.copy() 
        rho[1:] = rhon[1:] - V_max*dt/dx*(rhon[1:] -rhon[0:-1]) 
        rho[0] = boundary_0
    #---------------------------------------
    plt.plot(numpy.linspace(0,L,nx), rho, color='#003366', ls='--', lw=3)
    plt.ylim(0,250);
    plt.xlim(0,L);     
    return rho


#------------MAIN FUNCTION-----------    


##############################
####-----PART A -----------
############################
    
#Model Parameters
V_max = 80.;   #[km/hr]
L = 11.;       #[km]
rho_max = 250.; #[cars/km]      

#Initialize v --- values from exercise!
nx = 51;
rho0 = numpy.ones(nx)*10    #initial value rho(0,x)
rho0[10:20] = 50            #----" " ----    
boundary = 10.0             #boundary value rho(t,0)
x = numpy.linspace(0,L,nx)
    
print "The minimal velocity at time t = 0 is. %f m/s" % (1000./60**2*min(V(traffic(0, V_max, L, rho_max, rho0, boundary, nx), V_max, rho_max))) 

print "Average velocity at t=3 min is %f [m/s]" % (1000./60**2*numpy.average(V(traffic(3, V_max, L, rho_max, rho0, boundary, nx), V_max, rho_max)))

print "Minimum velocity at t=6 min is %f [m/s]" % (1000./60**2*min(V(traffic(6, V_max, L, rho_max, rho0, boundary, nx), V_max, rho_max)))    
    
    
    
##############################
####-----PART B -----------
############################    

#Model Parameters
V_max = 136.;   #[km/hr]
L = 11.;       #[km]
rho_max = 250.; #[cars/km]      

#Initialize v --- values from exercise!
nx = 51;
rho0 = numpy.ones(nx)*20    #initial value rho(0,x)
rho0[10:20] = 50            #----" " ----    
boundary = 20.0             #boundary value rho(t,0)
    
print "The minimal velocity at time t = 0 is. %f m/s" % (1000./60**2*min(V(traffic(0, V_max, L, rho_max, rho0, boundary, nx), V_max, rho_max)))   

print "Average velocity at t=3 min is %f [m/s]" % (1000./60**2*numpy.average(V(traffic(3, V_max, L, rho_max, rho0, boundary, nx), V_max, rho_max)))

print "Minimum velocity at t=3 min is %f [m/s]" % (1000./60**2*min(V(traffic(3, V_max, L, rho_max, rho0, boundary, nx), V_max, rho_max)))    
 
    
        