#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 16:34:57 2023

@author: bhaskar
"""

import numpy as np
import math



def curtis(r_1_vec, r_2_vec, delta_t,mu,senseOfMotion):
    
    #mu = 1.32712440018e20 * (86400**2) / 149597870700.0**3# GM_sun, au^3/day^2
    #calculate delta_theta
    r_1 = np.linalg.norm(r_1_vec, ord=2)
    r_2 = np.linalg.norm(r_2_vec, ord=2)
    r1_cross_r2 = np.cross(r_1_vec,r_2_vec)
    prograde = 0
    #for prograde
    if (senseOfMotion == "prograde"):
        if(r1_cross_r2[2] >= 0):
            delta_theta = math.acos(np.dot(r_1_vec,r_2_vec)/( r_1*r_2 ))
        else:
            delta_theta = 2*np.pi - math.acos(np.dot(r_1_vec,r_2_vec)/( r_1*r_2 ))
    
    #retrograde
    else:
        if(r1_cross_r2[2] < 0):
            delta_theta = math.acos(np.dot(r_1_vec,r_2_vec)/( r_1*r_2 ))
        else:
            delta_theta = 2*np.pi - math.acos(np.dot(r_1_vec,r_2_vec)/( r_1*r_2 ))
    
    #Calculate A
    A = math.sin(delta_theta)* math.sqrt( r_1*r_2/(1-math.cos(delta_theta)) )
    
    #Calculate z
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
    def y(z):
        u = r_1+r_2+ A*(z*S(z)-1)/math.sqrt( C(z) )
        return u
    def F(z):
        u = (y(z)/C(z))**1.5 *S(z) + A*math.sqrt(y(z)) - math.sqrt(mu)*delta_t
        return u
    def F_prime(z):
        if (z!=0):
            u = (y(z)/C(z))**1.5 * ( (1/(2*z))*(C(z) - 1.5*S(z)/C(z)) + (3/4)*S(z)**2/C(z) ) \
                + (A/8)*( (3*S(z)/C(z))*math.sqrt(y(z)) + A*math.sqrt(C(z)/y(z)) )
        else:
            u = (math.sqrt(2)/40)*y(0)**1.5 + (A/8)*(math.sqrt(y(0))+A/math.sqrt(2*y(0)) )
        return u
    
    #Find z
    z = 0.0
    #print("\nworking\n")
    zvalues = np.linspace(-1.0,1.0,10)
    fvalues = np.zeros(np.size(zvalues))
    for i,zval in enumerate(zvalues):
        #print("yes")
        fvalues[i] = F(zval)
        
    #print("working")
    for i in range(np.size(zvalues)-1):
        if (np.sign(F(zvalues[i]))!=np.sign(F(zvalues[i+1]))):
            #print("if pass")
            z = zvalues[i]
    
    #print(z)
    tolerance = 1.e-8
    ratio = 1.0
    while (tolerance<ratio):
        ratio = F(z)/F_prime(z)
        z = z - ratio
    #print (z)
    #Find Lagrange f and g functions
    f = 1 - y(z)/r_1
    g = A*math.sqrt(y(z)/mu)
    #f_dot = (math.sqrt(mu)/(r_1*r_2))* math.sqrt(y(z)/C(z))* (z*S(z)-1)
    g_dot = 1 - y(z)/r_2
    
    #Calculate v_1_vec and v_2_vec required for the trajectory of satellite
    v_1_vec = (1/g)*(r_2_vec-f*r_1_vec)
    v_2_vec = (1/g)*(g_dot*r_2_vec-r_1_vec)
    
    return [v_1_vec,v_2_vec]