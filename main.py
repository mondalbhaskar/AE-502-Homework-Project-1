#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 16:34:36 2023

@author: bhaskar
"""
#http://www.csgnetwork.com/julianmodifdateconv.html MJD converter
#import other files containing functions
import universal2BodyPropagator
import lambert

#import necessary modules
import numpy as np
import math
import matplotlib.pyplot as plt


#initial values provided
r_0_vec_earth = np.array([-1.796136509111975e-1, 9.667949206859814e-1, \
                    -3.668681017942158e-5]) #au
#r_0_vec = [1,2,3]
v_0_vec_earth = np.array([-1.720038360888334e-2, -3.211186197806460e-3, \
                    7.927736735960840e-7])#au/day

mu = 1.32712440018e20 * (86400**2) / 149597870700.0**3# GM_sun, au^3/day^2

target = 2 #Select target object here
mode = "rendezvous" # flyby/rendezvous

if (target==1):
    #1I/’Oumuamoua
    r_0_vec_target = np.array([3.515868886595499e-2, -3.162046390773074, \
                               4.493983111703389]) #au
    v_0_vec_target = np.array([-2.317577766980901e-3, 9.843360903693031e-3, \
                               -1.541856855538041e-2]) #au/day
    departureDates = np.arange(2457754.5,2458118.5,1.0) #Jan 1, 2017 - Dec 31, 2017
    arrivalDates = np.arange(2457966.5,2458514.5,1.0) #Aug 1, 2017 - Jan 31, 2019
    if (mode == "rendezvous"):
        max_delta_v = 50.0 # km/s
    else:
        max_delta_v = 20.0 # km/s
    
if (target==2):
    #2I/Borisov
    r_0_vec_target = np.array([7.249472033259724, 14.61063037906177, \
                               14.24274452216359]) #au
    v_0_vec_target = np.array([-8.241709369476881e-3, -1.156219024581502e-2, \
                               -1.317135977481448e-2]) #au/day
    
    departureDates = np.arange(2457754.5,2459061.5,1.0) #Jan 1, 2017 - Jul 31, 2020
    arrivalDates = np.arange(2458635.5,2459610.5,1.0) #Jun 1, 2019 - Jan 31, 2022
    if (mode == "rendezvous"):
        max_delta_v = 60.0 # km/s
    else:
        max_delta_v = 20.0 # km/s

T0 = 2457754.5 #Julian date for all the initial r_vec and v_vec values provided

#result arrays
#dayOfTakeoff_array=[]
#dayOfArrival_array=[]
#delta_t_array = []
dayOfArrival_array = np.empty([np.size(arrivalDates),np.size(departureDates)])
dayOfTakeoff_array = np.empty([np.size(arrivalDates),np.size(departureDates)])
delta_t_array=np.empty([np.size(arrivalDates),np.size(departureDates)])
#delta_v_array=np.empty([np.size(departureDates)*np.size(arrivalDates), np.size(departureDates)*np.size(arrivalDates)])
delta_v_array=np.empty([np.size(departureDates),np.size(arrivalDates)])
delta_v_array_flyby=np.empty([np.size(departureDates),np.size(arrivalDates)]) ##
cases=0
departs=0
print("Total dep: ", np.size(departureDates))
for i,departureDate in enumerate(departureDates):
    print("departures: ",i)
    departureTimeFromBase = departureDate - T0
    if (departureTimeFromBase==0):
        r_vec_start = r_0_vec_earth
        v_vec_start = v_0_vec_earth
    else:
        #use universal@BodyPropagator code
        [r_vec_start,v_vec_start] = universal2BodyPropagator.propagate(r_0_vec_earth,v_0_vec_earth,departureTimeFromBase,mu)
    
    for j,arrivalDate in enumerate(arrivalDates):
        arrivalTimeFromBase = arrivalDate - T0
        delta_t = arrivalDate-departureDate
        if (delta_t>0):
            [r_vec_end,v_vec_end] = universal2BodyPropagator.propagate(r_0_vec_target,v_0_vec_target,arrivalTimeFromBase,mu)
            
            try:
                [v_vec_start_reqd_prograde,v_vec_end_reqd_prograde] = lambert.curtis(r_vec_start,r_vec_end,delta_t,mu,"prograde")
            except:
                [v_vec_start_reqd_prograde,v_vec_end_reqd_prograde] = [[np.nan,np.nan,np.nan],[np.nan,np.nan,np.nan]]
            
            try:
                [v_vec_start_reqd_retrograde,v_vec_end_reqd_retrograde] = lambert.curtis(r_vec_start,r_vec_end,delta_t,mu,"prograde")
            except:
                [v_vec_start_reqd_retrograde,v_vec_end_reqd_retrograde] = [[np.nan,np.nan,np.nan],[np.nan,np.nan,np.nan]]
            
            delta_v_start_prograde = np.linalg.norm(v_vec_start_reqd_prograde - v_vec_start, ord=2)
            delta_v_start_retrograde = np.linalg.norm(v_vec_start_reqd_retrograde - v_vec_start, ord=2)
            if (mode=="rendezvous"):
                delta_v_end_prograde = np.linalg.norm(v_vec_end_reqd_prograde - v_vec_end, ord=2)
                delta_v_end_retrograde = np.linalg.norm(v_vec_end_reqd_retrograde - v_vec_end, ord=2)
            else:
                delta_v_end_prograde = 0
                delta_v_end_retrograde = 0
            
            delta_v_prograde = delta_v_start_prograde + delta_v_end_prograde
            delta_v_retrograde = delta_v_start_retrograde + delta_v_end_retrograde
            delta_v = (min(delta_v_prograde,delta_v_retrograde)) * 149597870.7/86400  #need to change this to km/s
            delta_v_flyby = (min(delta_v_start_prograde,delta_v_start_retrograde)) * 149597870.7/86400##
            if (delta_v < max_delta_v):
                delta_v_array[i,j]=delta_v
                #dayOfTakeoff_array[j,i] = departureDate -T0
                #dayOfArrival_array[j,i] = arrivalDate -arrivalDates[0]
                #delta_t_array[j,i] = delta_t
            else:
                delta_v_array[i,j]=np.nan
                #dayOfTakeoff_array[j,i] = departureDate -T0
                #dayOfArrival_array[j,i] = arrivalDate -arrivalDates[0]
                #delta_t_array[j,i] = delta_t
            if (delta_v_flyby<20.0):##
                delta_v_array_flyby[i,j]=delta_v_flyby
            else:
                delta_v_array_flyby[i,j]=np.nan
        else:
            delta_v_array[i,j]=np.nan # arrival date is before launch date
            delta_v_array_flyby[i,j]=np.nan ##
        dayOfTakeoff_array[j,i] = departureDate -T0
        dayOfArrival_array[j,i] = arrivalDate -arrivalDates[0]
        delta_t_array[j,i] = delta_t
    



Z = np.transpose(np.array(delta_v_array))
plt.figure(1)
plt.contourf(dayOfTakeoff_array, dayOfArrival_array, Z, levels=100,cmap='jet')

# Add labels and title
plt.xlabel('Day of launch (Days since Jan 1, 2017)')
plt.ylabel('Day of arrival (Days since Jun 1, 2019)')

plt.title('Contour Plot of $\Delta v$ for rendezvous with 2I/Borisov')
plt.colorbar(label='$\Delta v$ (km/s)') # Add a colorbar to a plot

# Show the plot
plt.savefig('2Irendezvous.png',bbox_inches='tight')




Z2 = np.transpose(np.array(delta_v_array_flyby))
plt.figure(2)
plt.contourf(dayOfTakeoff_array, dayOfArrival_array, Z2, levels=100,cmap='jet')

# Add labels and title
plt.xlabel('Day of launch (Days since Jan 1, 2017)')
plt.ylabel('Day of arrival (Days since Jun 1, 2019)')

plt.title('Contour Plot of $\Delta v$ for flyby of 2I/Borisov')
plt.colorbar(label='$\Delta v$ (km/s)') # Add a colorbar to a plot

# Show the plot
plt.savefig('2Iflyby.png',bbox_inches='tight')
plt.show()


