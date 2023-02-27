#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 01:11:58 2023

@author: bhaskar
"""

#import other files containing functions
#import universal2BodyPropagator
import lambert

#import necessary modules
import numpy as np
import math
import matplotlib.pyplot as plt

r1 = np.array([5644.0, 2830.0, 4170.0]) # km
r2 = np.array([-2240.0, 7320.0, 4980.0]) # km
mu = 3.986004418e5 #km^3/s^2
tof=1200 #s


[v_vec_start_reqd,v_vec_end_reqd] = lambert.curtis(r1,r2,tof,mu,"retrograde")

