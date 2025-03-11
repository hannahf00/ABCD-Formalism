# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 11:25:29 2025

@author: hanna
"""


'''
Calculates the waist size in the center of the crystal, if you give it the 
incoupling waist!


'''


import numpy as np

lc = 12         # length of crystal in mm
n = 1.66        # refr. index of crystal
l = 406         # total length 

l1 = 61.5      # length between the curved mirrors
l2 = l-l1      # length between all the remaining mirrors
l3 = (l-lc)/2   # length from incoupling position until first curved mirror
l4 = (l1-lc)/2  # length between first curved mirror and crystal surface
r = 50          # ROC of mirrors
lam = 902e-6    # wavelength
#w1 = 23.5e-3    # waist in crystal center
w2 = 206e-3     # incoupling waist
theta = 13.5/2*np.pi/180

cry = np.array([[1, lc/(2*n)], 
                [0, 1]])

short = np.array([[1, l4], 
                  [0, 1]])

mirr = np.array([[1, 0], 
                 [-2/(r*np.cos(theta)), 1]])

long = np.array([[1, l3], 
                 [0, 1]])

matrices = [cry, short, mirr, long
     
    ]

abcd = np.linalg.multi_dot(matrices)

print(abcd)
 
a = abcd[0, 0]
b = abcd[0, 1]
c = abcd[1, 0]
d = abcd[1, 1]

q2 = np.pi*w2**2/(-1j*lam)  # incoupling q-parameter

q1 = (a*q2+b)/(c*q2+d)      # q-paramter at crystal center

w1 = np.sqrt(-1j*lam*q1/np.pi)

print(np.abs(w1))
