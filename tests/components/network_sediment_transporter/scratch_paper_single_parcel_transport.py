# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 12:36:07 2019

@author: pfeif
"""

import numpy as np

rho = 1000
rho_s = 2650
S = 0.01
h = 2 
D = 0.05
parcel_vol = 1
link_length = 100 
dt = 60 # 1 min
g = 9.8
w = 15
R = (rho_s - rho)/rho

Tau = rho*g*h*S
Tau_star = Tau/((rho_s-rho)*g*D)

L_a= 0.515*D*(3.09*(Tau_star - 0.0549)**0.56)

capacity = L_a*link_length*w

##

Tau_star_ref = 0.021 + 0.015*np.exp(-20*0)
Tau_ref = Tau_star_ref * (rho_s - rho)*g*D

Tau_Taur_ratio = Tau / Tau_ref

W = 14*(1-(0.894/(Tau_Taur_ratio**0.5)))**4.5

p_vel = W*(Tau**(3/2))/((rho**(3/2))*g*R*L_a)

distance_to_travel = p_vel * dt

timestep = np.arange(0,12)

total_distance_traveled = distance_to_travel*timestep

final_location_in_link = (total_distance_traveled[-1]-200)/100