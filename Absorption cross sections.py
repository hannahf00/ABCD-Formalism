# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 09:35:30 2025

@author: hanna
"""

'''
ABSORPTION CROSS SECTIONS

This Program gives the answer to the question: How much power at which Waist 
size do we need to make diff. metal clusters absorb 5-6 photons. 

Metal clusters: Na, Hf, Au, Si 
For two different scenarios: 
    1) m = 100kDa and v = 130m/s
    2) m = 1MDa and v = 30m/s
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as c

wavelength = 226e-9 # in m
nu = c.c/wavelength

'Electric permittivity of diff. clusters'
'source: refractiveindex.info'
e_Na = -1.01 + 0.09899j     # Sodium, lam = 312nm
e_Hf = 0.516 + 3.3082j      # Hafnium, lam = 121nm 
e_Au = -0.39802 + 3.7963j   # Gold, lam = 225nm
e_Si = -9.0691 + 8.7624j    # Silicon, lam = 225nm


'Cluster masses'
m_1 =  100e3 * c.u   # 100 kDa
m_2 = 1e6 * c.u      # 1 MDa

'Cluster velocities'
v1 = 130 # in m/s ... for the 100kDa clusters
v2 = 30 # in m/s ... for the 1MDa clusters

'Densities in kg/m^3'
rho_Na = 9.7e2
rho_Hf = 1.3281e4 
rho_Au = 1.932e4
rho_Si = 2.33e3


def sigma_abs(e, rho, m1 = m_1, m2 = m_2, lam = wavelength): 
    ''' Calculates the absorption coefficient for diff. materials.
    
    lam:    wavelength
    e:      el. permittivity
    m:      mass in kg
    rho:    density in kg/m^3
    
    with r**3 = 3/4 * m/(np.pi*rho)
    
    RETURNS: absorption cross section in m^2 SI units for two different masses
    m1 = 100kDa and
    m2 = 1MDa
    
    '''
    
    sigma_abs_m1_SI = (8 * np.pi**2 * (3/4*m1/(np.pi*rho)) ) /lam * np.imag( (e-1)/(e+2) ) # in m^2 SI units
    #sigma_abs_m1 = sigma_abs_m1_SI * 1e4 # in cm^2
    
    sigma_abs_m2_SI = (8 * np.pi**2 * (3/4*m2/(np.pi*rho)) ) /lam * np.imag( (e-1)/(e+2) ) # in m^2 SI units
    #sigma_abs_m2 = sigma_abs_m2_SI * 1e4 # in cm^2
    
    return(sigma_abs_m1_SI,  sigma_abs_m2_SI)

#print(sigma_abs(e_Na, rho_Na))


' Calculated absorption cross sections in SI units '
sigma_abs_Na = sigma_abs(e_Na, rho_Na)
sigma_abs_Hf = sigma_abs(e_Hf, rho_Hf)
sigma_abs_Au = sigma_abs(e_Au, rho_Au)
sigma_abs_Si = sigma_abs(e_Si, rho_Si)


' Collecting the data for the 2 diff. scenarios '
Na_1 = np.array([sigma_abs_Na[0], v1])
Na_2 = np.array([sigma_abs_Na[1], v2])

Hf_1 = np.array([sigma_abs_Hf[0], v1])
Hf_2 = np.array([sigma_abs_Hf[1], v2])

Au_1 = np.array([sigma_abs_Au[0], v1])
Au_2 = np.array([sigma_abs_Au[1], v2])

Si_1 = np.array([sigma_abs_Si[0], v1])
Si_2 = np.array([sigma_abs_Si[1], v2])


def required_power(element, w_0, n = 10, nu = c.c/wavelength): 
    '''
    Calculates the relationship between Power and Waist size P(w) to absorb n 
    number of photons. The interaction time was calcualted over the mean 
    velocity and the waist size. The calculated power takes into account, that
    it is a standing wave by a factor of 4. Because the total electric field 
    is doubled, since the wave is traveling back and forth and the power is 
    proportional to the square of the el. field: P prop. (2*E)^2 
    --> this factor of 4 cancels in our case, bec. we have 3 depletion gratings
    and one beam for ionization 
    
    sigma: np.array that contains the absorption cross section for diff. masses
    n: nr. of photons to be absorbed 
    nu: frequency
    
    RETURNS: 
        The required power for different waist sizes. 
    '''
    
    sigma = element[0]  # extracting the corss section from the array
    v = element[1]      # extracting the corresponding velocity from the array
    
    P = n * c.h*nu/sigma * np.pi*w_0/np.sqrt(np.pi)* v  # power for velocity 1
    
    return(P)
   
    
def plot_m1_v1(): 
    w_0 = np.linspace(100e-6, 300e-6, 2) # waist from 100 um to 300 um 
    
    'FIGURE for smaller clustes: m = 100kDa and v = 130m/s'
    P_Na_1 = required_power(Na_1, w_0)
    P_Hf_1 = required_power(Hf_1, w_0)
    P_Au_1 = required_power(Au_1, w_0)
    P_Si_1 = required_power(Si_1, w_0)
    
    fig, ax = plt.subplots()
  
    ax.plot(w_0, P_Na_1, 'gold', label = 'Na') 
    ax.plot(w_0, P_Hf_1, 'dodgerblue', label = 'Hf')
    ax.plot(w_0, P_Au_1, 'forestgreen', label = 'Au')
    ax.plot(w_0, P_Si_1, 'firebrick', label = 'Si')
    
    ax.set_xticks(np.linspace(100e-6, 300e-6, 5))  # Adjust ticks
    ax.set_xticklabels([f"{tick * 1e6:.0f}" for tick in ax.get_xticks()])
    
    ax.set_yticks(np.linspace(0, 1.4, 8))  # Adjust ticks
    ax.set_yticklabels([f"{tick * 1e3:.0f}" for tick in ax.get_yticks()])
    
    plt.xlabel('waist size ($\mu$m)')
    plt.ylabel('laser power (mW)')
    plt.title('Required Power for Photon Absorption \n mass = 100kDa, v = 130m/s')
    ax.legend()
    ax.grid(True)
    plt.rcParams['figure.dpi']=2000
    plt.show
   
def plot_m2_v2(): 
    
    'FIGURE for bigger clustes: m = 1MDa and v = 30m/s'
    
    w_0 = np.linspace(100e-6, 300e-6, 2) # waist from 100 um to 300 um 
    
    P_Na_2 = required_power(Na_2, w_0)
    P_Hf_2 = required_power(Hf_2, w_0)
    P_Au_2 = required_power(Au_2, w_0)
    P_Si_2 = required_power(Si_2, w_0)
    
     
    fig, ax1 = plt.subplots()
  
    ax1.plot(w_0, P_Na_2, 'gold', label = 'Na') 
    ax1.plot(w_0, P_Hf_2, 'dodgerblue', label = 'Hf')
    ax1.plot(w_0, P_Au_2, 'forestgreen', label = 'Au')
    ax1.plot(w_0, P_Si_2, 'firebrick', label = 'Si')
    
    ax1.set_xticks(np.linspace(100e-6, 300e-6, 5))  # Adjust ticks
    ax1.set_xticklabels([f"{tick * 1e6:.0f}" for tick in ax1.get_xticks()])
    
    ax1.set_yticks(np.linspace(0, 0.027, 9))  # Adjust ticks
    ax1.set_yticklabels([f"{tick * 1e3:.0f}" for tick in ax1.get_yticks()])
    
    plt.xlabel('waist size ($\mu$m)')
    plt.ylabel('laser power (mW)')
    plt.title('Required Power for Photon Absorption \n mass = 1MDa, v = 30m/s')
    ax1.legend()
    ax1.grid(True)
    plt.rcParams['figure.dpi']=2000
    plt.show
    
    
    
    #plt.show
    
plot_m1_v1()
#plot_m2_v2()

    


    

