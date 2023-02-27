#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 12:58:50 2023

@author: bhaskar
"""

#import other files containing functions
import universal2BodyPropagator
#import lambert

#import necessary modules
import numpy as np
import math
import matplotlib.pyplot as plt

mu = 1.32712440018e20 * (86400**2) / 149597870700.0**3# GM_sun, au^3/day^2

#initial values provided
r_0_vec = np.array([-1.793111147418941E-01, 9.685804332587541E-01, \
                    1.612991618128246E-04]) #au
#r_0_vec = [1,2,3]
v_0_vec = np.array([-1.721950145453741E-02, -3.058649113379992E-03, \
                    -3.172877676994421E-07])#au/day

delta_t = 120.0

[r_1_vec,v_1_vec] = universal2BodyPropagator.propagate(r_0_vec,v_0_vec,delta_t,mu)
