# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 18:43:38 2023

@author: Mika
"""
import math


### -------------------------------------------------------------------
# =============================================================================
# # Question 1
# =============================================================================

# Parameters:
a = 1.16/6.53
p_ice = 917
p_water = 1000
E=8.7*10**9
h_ice = 0.28
g= 9.81
mu = 0.33
r = 1.16
rho_f =1.5*10**6
x = 1
y = 4.25
n = 1- p_ice/p_water

## consideration 1: -----------------
mass_1 = (p_water - p_ice)*h_ice*x*y
print(f'mass for assumption 1 = {round(mass_1, 2)} kg')

## consideration 2: -----------------

Lambda = ((E*h_ice**3)/(12*p_water*g*(1-mu**2)))**(1/4)
    
c0 =-1/a+1/8*math.pi*a-1/16*(1.3659-math.log(a))*a**3
Pw = (n*h_ice*p_water*g)/(1+a*c0)*math.pi*r**2

mass_2 = Pw/g
print(f'mass for assumption 2 = {round(mass_2)} kg')

## consideration 3:
c1 = ((0.619-math.log(a))*a/2+math.pi*a**3/64)/(math.pi*a)
P_cr = (rho_f*h_ice**2)/(3*(1+mu)*c1)
P_f_u = 2*P_cr
mass_3=P_f_u/g
print(f'mass for assumption 3 = {round(mass_3)} kg')

    
# =============================================================================
# # Question 2
# =============================================================================
alpha = 5.0*10**-5
distance_b = 1.42*10**3 #m
distance_a = 3*10**3 #m
degreeincrease = 7
extension_a = alpha*distance_a*degreeincrease
print(f'extension in direction a: {round(extension_a, 2)} m')
extension_b = alpha*distance_b*degreeincrease
print(f'extension in direction b: {extension_b} m')
