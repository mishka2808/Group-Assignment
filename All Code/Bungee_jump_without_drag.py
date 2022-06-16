"""
Created on Tues Jun 14 11:47:23 2022

@author: Mish & Dhers
//Graphing from the given data (Without Drag)//
1. Beta Version - Mish
2. Final Version - Dhers
"""
import pandas as pd
import os
import numpy as np

path = os.getcwd()
name = '/zero_drag'

df = pd.read_csv(path+name+'.csv')      #loads the called .csv file from the local downloads folder

x = np.asarray(df.x)        #takes the x positions from what was loaded, puts them into an array

y = np.asarray(df.y)        #puts y positions into array

vx = np.asarray(df.vx)      #puts x velocities into array

vy = np.asarray(df.vy)      #puts y velocities into array


Ekin = np.asarray(df.E_kin)      #puts kinetic energies into array

Epot = np.asarray(df.E_pot)      #puts potential energies into array

Etot = np.asarray(df.E_total)    #puts total energies into array

import matplotlib.pyplot as plt

plt.figure()
plt.plot([-20], [0], "or", label='A', markersize=10)
plt.plot([20], [0], "ob", label='B',markersize=10)
plt.plot(x,y,"ok")
plt.xlabel('X position')
plt.ylabel('Y position')
plt.legend()
plt.show()
#plot of the path travelled by the object

plt.figure()
plt.plot(Etot, '-k', label = 'Total Energy')
plt.plot(Epot, '-b', label = 'Potential Energy')
plt.plot(Ekin, '-r', label = 'Kinetic Energy')
plt.xlabel('Steps')
plt.ylabel('Energy')
plt.legend()
plt.show()
#plots the energies of the object as it travels along the path