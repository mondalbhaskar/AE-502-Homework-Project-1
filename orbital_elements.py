#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 16:10:53 2023

@author: bhaskar
"""
#import necessary modules
import numpy as np
import math
import matplotlib.pyplot as plt


mu = 1.32712440018e20 * (86400**2) / 149597870700.0**3# GM_sun, au^3/day^2

target = 1 #Select target object here

if (target==1):
    #1I/’Oumouamoua
    r_vec = np.array([3.515868886595499e-2, -3.162046390773074, \
                               4.493983111703389]) #au
    v_vec = np.array([-2.317577766980901e-3, 9.843360903693031e-3, \
                               -1.541856855538041e-2]) #au/day
  
if (target==2):
    #2I/Borisov
    r_vec = np.array([7.249472033259724, 14.61063037906177, \
                               14.24274452216359]) #au
    v_vec = np.array([-8.241709369476881e-3, -1.156219024581502e-2, \
                               -1.317135977481448e-2]) #au/day
    
r = np.linalg.norm(r_vec,ord=2)
v = np.linalg.norm(v_vec,ord=2)
h_vec = np.cross(r_vec,v_vec)
h = np.linalg.norm(h_vec,ord=2) # angular momentum

e_vec = np.cross(v_vec,h_vec)/mu - r_vec/r
e = np.linalg.norm(e_vec,ord=2) #eccentricity
print('eccentricity = ',e)