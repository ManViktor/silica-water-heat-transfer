# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 15:02:40 2023

@author: mandrolk1
"""

from pylab import *
import numpy as np
import matplotlib.pyplot as plt


def angle_between_points(P1, P2, P3):
    vector1 = (P1[0] - P2[0], P1[1] - P2[1], P1[2] - P2[2])
    vector2 = (P3[0] - P2[0], P3[1] - P2[1], P3[2] - P2[2])
    dot_product = vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]
    magnitude1 = math.sqrt(vector1[0] ** 2 + vector1[1] ** 2 + vector1[2] ** 2)
    magnitude2 = math.sqrt(vector2[0] ** 2 + vector2[1] ** 2 + vector2[2] ** 2)
    angle_rad = math.acos(dot_product / (magnitude1 * magnitude2))
    angle_deg = math.degrees(angle_rad)
    return angle_deg

def distance_between_points(P1, P2):
    vector = (P1[0] - P2[0], P1[1] - P2[1])
    magnitude = math.sqrt(vector1[0] ** 2 + vector1[1] ** 2)
    return magnitude





data = loadtxt('C:/LEMTA/Data/SiO2/ClayFF/creation/SiO2_alpha_orthogonal_cell.data', skiprows=17)

x_cell = 5.027782000000002
y_cell = 8.708374
z_cell = 5.518918000000006



x = data[:,4]
y = data[:,5]
z = data[:,6]
atype = data[:,2]
charge = data[:,3]


x -= min(x)      #shift
x_cell -= min(x)



mult_x = 12
mult_y = 16
mult_z = 10
X,Y,Z, new_atype, new_charge = [], [], [], [], []
for i in range(mult_x):
    for j in range(mult_y):
        for k in range(mult_z):
            for a in range(len(x)):
                if z[a]+z_cell*k < 54:
                    X.append(x[a]+x_cell*i)
                    Y.append(y[a]+y_cell*j)
                    Z.append(z[a]+z_cell*k)
                    new_atype.append(atype[a])
                    if atype[a] == 1.0:
                        new_charge.append(2.1)
                    elif atype[a] == 2.0:
                        new_charge.append(-1.05)
                

Si_x = array(X)
Si_y = array(Y)
Si_z = array(Z)

Si = []
for i in range(len(Si_x)):
    Si.append([Si_x[i], Si_y[i], Si_z[i]])
Si = array(Si)    




 
    
#Writing
output_file = 'C:/LEMTA/Data/SiO2/ClayFF/creation/system_configurations/ITR/SiO2_base.data'
with open(output_file, 'w') as fdata:
    fdata.write('##SiO2 LAMMPS slab\n\n')
    fdata.write('{} atoms\n\n'.format(len(X)))
    
    fdata.write('{} atom types\n\n'.format(2))
    
    fdata.write('{} {} xlo xhi\n'.format(0.0, max(X)+x_cell-max(x)))
    fdata.write('{} {} ylo yhi\n'.format(0.0, max(Y)+y_cell-max(y)))
    fdata.write('{} {} zlo zhi\n'.format(0.0, max(Z)+z_cell-max(z))) 
    fdata.write('\n')
    
    fdata.write('Masses\n\n')
    fdata.write('1 28.0855 # Si\n')
    fdata.write('2 15.9990 # O\n\n')

    
    
    
    #Atoms section
    fdata.write('Atoms #full\n\n')

    for i, pos in enumerate(X):
        fdata.write('{} 1 {} {} {} {} {}\n'.format(i+1, int(new_atype[i]), new_charge[i], Si_x[i], Si_y[i], Si_z[i] ))
    fdata.close()