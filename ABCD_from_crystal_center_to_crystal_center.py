# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 09:42:28 2025

@author: hanna
"""

'''
This is just a test for self-consistency! It calculates the ABCD-round trip 
matrix for our IR-VIS Cavity. 

Starting point = end point = center of crystal

It shows that the q-paramter stays the same after one full round trip!

'''




import numpy as np

lc = 12         # length of crystal in mm
n = 1.66        # refr. index of crystal
l1 = 60.8       # length between the curved mirrors
l2 = 345.2      # length between all the remaining mirrors
r = 50          # ROC of mirrors
lam = 902e-6    # wavelength
w1 = 23.5e-3    # waist in crystal center

cry = np.array([[1, lc/(2*n)], 
                [0, 1]])

short = np.array([[1, (l1-lc)/2], 
                  [0, 1]])

mirr = np.array([[1, 0], 
                 [-2/r, 1]])

long = np.array([[1, l2], 
                 [0, 1]])

matrices = [
     cry, short, mirr, long, mirr, short, cry
    ]

abcd = np.linalg.multi_dot(matrices)


print(abcd)
 
a = abcd[0, 0]
b = abcd[0, 1]
c = abcd[1, 0]
d = abcd[1, 1]



q1 = np.pi*w1**2/(-1j*lam)

q2 = (a*q1+b)/(c*q1+d)

print(q1, q2)
