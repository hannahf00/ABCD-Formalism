# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 18:03:19 2025

@author: hanna
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameter
lambda_ = 902e-9    # wavelength in m
n_air = 1           # refr. index of air
n_crystal = 1.66    # refr. idex of crystal
R = 50e-3           # ROC of focusing mirrors 50 mm in m
theta = 13.5/2      # in degree

# waist at center of crystal
w0_h = 24.3e-6  # horizontal 
w0_v = 25.7e-6  # vertical 

# Rayleigh-Range in crystal 
zR_h = np.pi * n_crystal * w0_h**2 / lambda_
zR_v = np.pi * n_crystal * w0_v**2 / lambda_

# q-Parameter at crystal center
q_h = 1j * zR_h
q_v = 1j * zR_v

# effective radius of curvature of mirrors, bec. of tilt
R_eff_h = R * np.cos(np.radians(theta))  # horizontal 
R_eff_v = R / np.cos(np.radians(theta))  # vertical 

# different segments
segments = [
    {"type": "propagation", "length": 6e-3, "n": n_crystal},
    {"type": "interface", "factor": 1/n_crystal, "n": n_air},
    {"type": "propagation", "length": 24.4e-3, "n": n_air},
    {"type": "mirror"},
    {"type": "propagation", "length": 345.2e-3, "n": n_air},
    {"type": "mirror"},
    {"type": "propagation", "length": 24.4e-3, "n": n_air},
    {"type": "interface", "factor": n_crystal, "n": n_crystal},
    {"type": "propagation", "length": 6e-3, "n": n_crystal}
]

dz = 1e-4  # stepsize
# arrays for storage 
z_all = []
w_h_all = []
w_v_all = []
z_current = 0

# calculates the waist propagation 
def calc_w(q, n_local):
    return np.sqrt(lambda_ / (np.pi * n_local * abs(np.imag(1/q))))

# propagation along one segment
for S in segments:
    if S["type"] == "propagation":
        L, n_local = S["length"], S["n"]
        d_vals = np.arange(0, L + dz, dz)
        for d in d_vals:
            q_h_sample, q_v_sample = q_h + d, q_v + d
            z_all.append(z_current + d)
            w_h_all.append(calc_w(q_h_sample, n_local))
            w_v_all.append(calc_w(q_v_sample, n_local))
        q_h += L
        q_v += L
        z_current += L
    
    elif S["type"] == "interface":
        factor = S["factor"]
        q_h *= factor
        q_v *= factor
        z_all.append(z_current)
        w_h_all.append(calc_w(q_h, S["n"]))
        w_v_all.append(calc_w(q_v, S["n"]))
    
    elif S["type"] == "mirror":
        q_h = q_h / (1 - (2/R_eff_h) * q_h)
        q_v = q_v / (1 - (2/R_eff_v) * q_v)
        z_all.append(z_current)
        w_h_all.append(calc_w(q_h, n_air))
        w_v_all.append(calc_w(q_v, n_air))

# converts lists to arrays
z_all = np.array(z_all) * 1e3
w_h_all = np.array(w_h_all) * 1e6
w_v_all = np.array(w_v_all) * 1e6

# plots the results
plt.figure(dpi=2000)
plt.plot(z_all, w_h_all, 'b-', linewidth=1.5, label='horizontal')
plt.plot(z_all, w_v_all, 'r-', linewidth=1.5, label='vertical')
plt.xlabel('distance z [mm]')
plt.ylabel('waist w(z) [µm]')
plt.legend(loc='best')
plt.title('Waist size w(z) for a full round-trip')
plt.grid(True)
plt.show()

