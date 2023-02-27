#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 15:56:58 2023

@author: bhaskar
"""

import numpy as np
import math

def f(r_0, v_r0, mu, chi, alpha,delta_t):
    u = r_0*v_r0* mu**(-0.5) * chi**2 * C(alpha*chi**2) + (1-alpha*r_0)* chi**3 \
        * S(alpha*chi**2) + r_0*chi - mu**(0.5) * delta_t
    return u
def f_prime(r_0, v_r0, mu, chi, alpha):
    u = r_0*v_r0* mu**(-0.5) * chi * (1-alpha*chi**2* S(alpha*chi**2)) \
        + (1-alpha*r_0)*chi**2*C(alpha*chi**2) + r_0
    return u


def C(z):
    if z > 0:
        u = (1-math.cos(math.sqrt(z))) / z 
    elif z < 0:
        u = (math.cosh(math.sqrt(-z))-1)/(-z)
    else:
        u = 0.5
    return u

def S(z):
    if z > 0:
        u = (math.sqrt(z)-math.sin(math.sqrt(z)))/(z**1.5) 
    elif z < 0:
        u = (math.sinh(math.sqrt(-z))-math.sqrt(-z))/((-z)**1.5)
    else:
        u = 1/6
    return u

#0. Define constants
#r_0_vec = np.array([-1.796136509111975e-1, 9.667949206859814e-1, \
#                    -3.668681017942158e-5]) #au
#r_0_vec = [1,2,3]
#v_0_vec = np.array([-1.720038360888334e-2, -3.211186197806460e-3, \
#                    7.927736735960840e-7])#au/day
#mu = 1.32712440018e20 * (86400**2) / 149597870700.0**3# GM_sun, au^3/day^2, 86164.09054

#delta_t = 365.25 # days
def propagate(r_0_vec,v_0_vec, delta_t,mu):
    r_0 = np.linalg.norm(r_0_vec, ord=2)
    v_0 = np.linalg.norm(v_0_vec, ord=2)
    v_r0 = np.dot(r_0_vec,v_0_vec)/r_0
    alpha = 2.0/r_0 - v_0**2 / mu
    
    #1. get initial estimate of chi
    chi = math.sqrt(mu)*abs(alpha)*delta_t #initial estimate
    tolerance = 1.e-8
    ratio = 1.0
    #print(chi)
    while ratio > tolerance:
        # 2. & 3. Calculate f/f'
        ratio = f(r_0, v_r0, mu, chi, alpha,delta_t)/f_prime(r_0, v_r0, mu, chi, alpha)
        #4. Update chi, when ratio is larger
        chi = chi - ratio
        #print(chi)
    #.5 Loop ends when ratio is less than tolerance
    
    #Get f and g values
    
    f_val = 1-(chi**2/r_0)*C(alpha*chi**2)
    g_val = delta_t - mu**(-0.5)*chi**3*S(alpha*chi**2)
    
    r_vec = f_val*r_0_vec + g_val*v_0_vec
    r = np.linalg.norm(r_vec,ord=2)
    
    f_dot_val = (mu**0.5 /(r*r_0))*(alpha*chi**3*S(alpha*chi**2)-chi)
    g_dot_val = 1-(chi**2/r)*C(alpha*chi**2)
    
    v_vec = f_dot_val*r_0_vec + g_dot_val*v_0_vec
    
    return [r_vec,v_vec]