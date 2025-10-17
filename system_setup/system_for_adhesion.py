# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 13:53:27 2023

@author: mandrolk1
"""



from pylab import *
import numpy as np
import re
from scipy.optimize import fsolve


output_file = 'C:/LEMTA/Data/SiO2/ClayFF/creation/system_configurations/adhesion/System_Configuration.data'
si_input = 'C:/LEMTA/Data/SiO2/ClayFF/creation/system_configurations/adhesion/SiO2_base.data'
water_input = 'C:/LEMTA/Data/SiO2/ClayFF/creation/system_configurations/adhesion/water.data'

x_pbc_dist = 0.34937302950
y_pbc_dist = 0.69989201838


with open(si_input, 'r') as f:
    lines = f.readlines()
    
def row_num(word):
    with open(si_input, 'r') as f:
        lines = f.readlines()
    
    for n, i in enumerate(lines):
        if re.search(r'\b' + word + r'\b', i):
            row = n 
            break
    return row

def row_num_H2O(word):
    with open(water_input, 'r') as f:
        lines = f.readlines()
    
    for n, i in enumerate(lines):
        if re.search(r'\b' + word + r'\b', i):
            row = n 
            break
    return row
    

def angle_between_points(P1, P2, P3):
    vector1 = (P1[0] - P2[0], P1[1] - P2[1], P1[2] - P2[2])
    vector2 = (P3[0] - P2[0], P3[1] - P2[1], P3[2] - P2[2])
    dot_product = vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]
    magnitude1 = math.sqrt(vector1[0] ** 2 + vector1[1] ** 2 + vector1[2] ** 2)
    magnitude2 = math.sqrt(vector2[0] ** 2 + vector2[1] ** 2 + vector2[2] ** 2)
    angle_rad = math.acos(dot_product / (magnitude1 * magnitude2))
    angle_deg = math.degrees(angle_rad)
    return angle_deg


def exchange_elements(array1, array2, num_elements):
    # Get random indices from both arrays
    indices1 = np.random.choice(len(array1), num_elements, replace=False)
    indices2 = np.random.choice(len(array2), num_elements, replace=False)

    # Exchange elements between the arrays
    array1[indices1], array2[indices2] = array2[indices2].copy(), array1[indices1].copy()


###############################################----SILICA----#################################################
ra = 3 #row number where is string 'N atoms'

first_atom_raw = row_num('Atoms')+3   
amount_of_atoms = [int(s) for s in re.findall(r'\b\d+\b', lines[ra-1])][0]


#atoms
data_atoms_Si = loadtxt(si_input, skiprows = first_atom_raw-1, max_rows=amount_of_atoms)

Si_x = data_atoms_Si[:,4]
Si_y = data_atoms_Si[:,5]
Si_z = data_atoms_Si[:,6]
Si_type = data_atoms_Si[:,2]
Si_charge = data_atoms_Si[:,3]
Si_id = data_atoms_Si[:,0]





Si, Atoms = [], []
for i in range(len(Si_x)):
    Si.append([Si_x[i], Si_y[i], Si_z[i], Si_id[i]])
    Atoms.append([Si_id[i], Si_type[i], Si_charge[i], Si_x[i], Si_y[i], Si_z[i]])
Si = array(Si)  

# create arrays with alterating upper Si atoms 
Si_up = Si[Si_z==unique(Si_z)[-1]] 
Si_down = Si[Si_z==unique(Si_z)[0]] 

np.random.shuffle(Si_up)

# Assuming 'Si_up' is your input array
# Adjust the percentage as needed
percentage_Si_up1 = 37.5  # Percentage for Si_up1
percentage_Si_up2 = 100.0 - percentage_Si_up1  # Percentage for Si_up2

# Determine the grid size and dimensions
grid_spacing = 1.5  # Adjust as needed
x_max = np.max(Si_up[:, 0])
x_min = np.min(Si_up[:, 0])
y_max = np.max(Si_up[:, 1])
y_min = np.min(Si_up[:, 1])
x_range = int((x_max - x_min) / grid_spacing) + 1
y_range = int((y_max - y_min) / grid_spacing) + 1

# Calculate the number of atoms in 'Si_up'
num_atoms = len(Si_up)

# Calculate the number of atoms for 'Si_up1' and 'Si_up2' based on the specified percentage
num_atoms_Si_up1 = int(num_atoms * (percentage_Si_up1 / 100))
num_atoms_Si_up2 = num_atoms - num_atoms_Si_up1

# Initialize arrays for Si_up1 and Si_up2
Si_up1 = [] # connected to (OH)2
Si_up2 = [] # connected to (CH3)2

print('Start: '+str(sum(array(Atoms)[:,2])))

# Initialize counters for the number of atoms in each group
count_Si_up1 = 0
count_Si_up2 = 0

# Loop through the atoms and assign them to Si_up1 or Si_up2 based on grid position
for atom in Si_up:
    x, y, _, _ = atom
    grid_x = int((x - x_min) / grid_spacing)
    grid_y = int((y - y_min) / grid_spacing)
    if (grid_x + grid_y) % 2 == 0 and count_Si_up1 < num_atoms_Si_up1:
        Si_up1.append(atom)
        count_Si_up1 += 1
    else:
        Si_up2.append(atom)
        count_Si_up2 += 1

# If the loop didn't assign enough atoms to Si_up1, fill the rest from Si_up2
while count_Si_up1 < num_atoms_Si_up1:
    Si_up1.append(Si_up2.pop(0))
    count_Si_up1 += 1


Si_up1 = array(Si_up1)    # connected to (OH)2
Si_up2 = array(Si_up2)    # connected to (CH3)2



# Change type and charge for Si (Si_up2) that contact CH3
for i in range(len(Atoms)):
    if Atoms[i][0] in [row[3] for row in Si_up2]:
        Atoms[i][1] = int(5)
        Atoms[i][2] = 1.05


print('Change Si in CH3: '+str(sum(array(Atoms)[:,2])))

#### add additional O atom below to make neutral
O1down_shift_from_Si = [-0.4827, 1.3789999, -0.713779] #got from xyz
O2down_shift_from_Si = [0.13341, -1.1679, -1.125859] #got from xyz
O_down = []
for i in range(len(Si_down)):
    O1x = Si_down[:,0][i]+O1down_shift_from_Si[0]
    O1y = Si_down[:,1][i]+O1down_shift_from_Si[1]
    O1z = Si_down[:,2][i]+O1down_shift_from_Si[2]
    O1id = int(len(Atoms)+i+1)
    O_down.append([O1id, 2, -1.05, O1x, O1y, O1z])


Atoms = Atoms + O_down


print('Add one O below: '+str(sum(array(Atoms)[:,2])))


##############################################Create additional atoms###################################
fake_bl_SiO = 1.62    # bond length
bl_OH = 1.0
bl_SiC = 1.85
bl_CH = 1.09
ang_SiOH = 109.47
ang_SiCH = 109.5
ang_HCH = 107.8
ang_CSiC = 110.0
O1, H1, O2, H2, C = [], [], [], [], []



#---------------------------------------- OH -----------------------------------
Bonds, Angles = [], []

O1_shift_from_Si = [0.945, -0.7, 1.1258] #got from xyz
H1_shift_from_O1  = [-0.29298808, 0.2929880, 0.9101186]
for i in range(len(Si_up1)):
    Ox1 = Si_up1[:,0][i]+O1_shift_from_Si[0]
    if Ox1 > max(Si_x)+x_pbc_dist:
        Ox1 -= max(Si_x)+x_pbc_dist
    if Si_up1[:,1][i]+O1_shift_from_Si[1] < 0:
        Oy1 = Si_up1[:,1][i]+O1_shift_from_Si[1]+max(Si_y)+y_pbc_dist
    else:
        Oy1 = Si_up1[:,1][i]+O1_shift_from_Si[1]
    Oz1 = Si_up1[:,2][i]+O1_shift_from_Si[2]
    Oid1 = int(len(Atoms))+i+1
    O1.append([Oid1, 3, -0.95,  Ox1, Oy1, Oz1])
    
    
    Hx1 = Ox1+H1_shift_from_O1[0]
    if Hx1 > max(Si_x)+x_pbc_dist:
        Hx1 -= max(Si_x)+x_pbc_dist
    Hy1 = Oy1+H1_shift_from_O1[1]
    Hz1 = Oz1+H1_shift_from_O1[2]
    Hid1 = int(len(Atoms))+i+1+len(Si_up1)
    
    H1.append([Hid1, 4, 0.425, Hx1, Hy1, Hz1])
    
    
    Bonds.append([i+1, 1, Oid1, Hid1])
    Angles.append([i+1, 1, Si_up1[:,3][i], Oid1, Hid1])


O2_shift_from_Si = [-0.952, 1.1073, 0.7137] #got from xyz
H2_shift_from_O2  = [0.303296, 0.14277250, 0.9421397] 
for i in range(len(Si_up1)):   

    Ox2 = Si_up1[:,0][i]+O2_shift_from_Si[0]
    Oy2 = Si_up1[:,1][i]+O2_shift_from_Si[1]
    Oz2 = Si_up1[:,2][i]+O2_shift_from_Si[2]
    Oid2 = int(len(Atoms))+i+1+2*len(Si_up1)
    O2.append([Oid2, 3, -0.95,  Ox2, Oy2, Oz2])
    
    
    Hx2 = Ox2+H2_shift_from_O2[0]
    Hy2 = Oy2+H2_shift_from_O2[1]
    Hz2 = Oz2+H2_shift_from_O2[2]
    Hid2 = int(len(Atoms))+i+1+3*len(Si_up1)
    
    
    H2.append([Hid2, 4, 0.425, Hx2, Hy2, Hz2])

    Bonds.append([len(Si_up1)+i+1, 1, Oid2, Hid2])
    Angles.append([len(Si_up1)+i+1, 1, Si_up1[:,3][i], Oid2, Hid2])
    
    
Atoms = Atoms + O1 + H1 + O2 + H2
print('Add OH: '+str(sum(array(Atoms)[:,2])))

#---------------------------------------- (CH3)2 ---------------------   

C1_shift_from_Si = [1.14268464, 0.9954349 , 1.06107548]
H11_shift_from_Si = [0.67159978, 1.9456448 , 1.312631]        #H11 - first H in the first CH3 group
H12_shift_from_Si = [2.07075966, 1.18156702, 0.52058348]
H13_shift_from_Si = [1.36170325, 0.4428307 , 1.97472718]

C2_shift_from_Si = [-1.14268464, -0.9954349 , 1.06107548]
H21_shift_from_Si = [-0.67159978, -1.9456448 , 1.312631]        
H22_shift_from_Si = [-2.07075966, -1.18156702, 0.52058348]
H23_shift_from_Si = [-1.36170325, -0.4428307 , 1.97472718]


C1, H11, H12, H13 = [], [], [], []
C2, H21, H22, H23 = [], [], [], []
for i in range(len(Si_up2)):   
    
    # first CH3
    
    Cx = Si_up2[:,0][i]+C1_shift_from_Si[0]
    if Cx > max(Si_x)+x_pbc_dist:
        Cx -= max(Si_x)+x_pbc_dist
    Cy = Si_up2[:,1][i]+C1_shift_from_Si[1]
    if Cy < 0:
        Cy += max(Si_y)+y_pbc_dist
    Cz = Si_up2[:,2][i]+C1_shift_from_Si[2]
    Cid = int(len(Atoms))+i+1
    C1.append([Cid, 6, -0.18,  Cx, Cy, Cz])
    
    H1x = Si_up2[:,0][i]+H11_shift_from_Si[0]
    if H1x > max(Si_x)+x_pbc_dist:
        H1x -= max(Si_x)+x_pbc_dist
    H1y = Si_up2[:,1][i]+H11_shift_from_Si[1]
    if H1y < 0:
        H1y += max(Si_y)+y_pbc_dist
    H1z = Si_up2[:,2][i]+H11_shift_from_Si[2]
    H1id = int(len(Atoms))+i+1+len(Si_up2)
    H11.append([H1id, 7, 0.06, H1x, H1y, H1z])
    
    
    H2x = Si_up2[:,0][i]+H12_shift_from_Si[0]
    if H2x > max(Si_x)+x_pbc_dist:
        H2x -= max(Si_x)+x_pbc_dist
    H2y = Si_up2[:,1][i]+H12_shift_from_Si[1]
    if H2y < 0:
        H2y += max(Si_y)+y_pbc_dist
    H2z = Si_up2[:,2][i]+H12_shift_from_Si[2]
    H2id = int(len(Atoms))+i+1+2*len(Si_up2)
    H12.append([H2id, 7, 0.06, H2x, H2y, H2z])
    
    
    H3x = Si_up2[:,0][i]+H13_shift_from_Si[0]
    if H3x > max(Si_x)+x_pbc_dist:
        H3x -= max(Si_x)+x_pbc_dist
    H3y = Si_up2[:,1][i]+H13_shift_from_Si[1]
    if H3y < 0:
        H3y += max(Si_y)+y_pbc_dist
    H3z = Si_up2[:,2][i]+H13_shift_from_Si[2]
    H3id = int(len(Atoms))+i+1+3*len(Si_up2)
    H13.append([H3id, 7, 0.06, H3x, H3y, H3z])
    

    Bonds.append([1000, 2, int(Si_up2[:,3][i]), Cid])    # Si - C Numbering will do later
    Bonds.append([1000, 3, Cid, H1id])    # C - H1
    Bonds.append([1000, 3, Cid, H2id])    # C - H2
    Bonds.append([1000, 3, Cid, H3id])    # C - H3
    
    
    Angles.append([1000, 2, Si_up2[:,3][i], Cid, H1id])    #Si - C - H1
    Angles.append([1000, 2, Si_up2[:,3][i], Cid, H2id])    #Si - C - H2
    Angles.append([1000, 2, Si_up2[:,3][i], Cid, H3id])    #Si - C - H3
    Angles.append([1000, 3, H1id, Cid, H2id])    #H1 - C - H2
    Angles.append([1000, 3, H1id, Cid, H3id])    #H1 - C - H3
    Angles.append([1000, 3, H2id, Cid, H3id])    #H2 - C - H3
    
    
    # second CH3

    C2x = Si_up2[:,0][i]+C2_shift_from_Si[0]
    C2y = Si_up2[:,1][i]+C2_shift_from_Si[1]
    if C2y < 0:
        C2y += max(Si_y)+y_pbc_dist
    C2z = Si_up2[:,2][i]+C2_shift_from_Si[2]
    C2id = int(len(Atoms))+i+1+4*len(Si_up2)
    C2.append([C2id, 6, -0.18,  C2x, C2y, C2z])
    
    H1x2 = Si_up2[:,0][i]+H21_shift_from_Si[0]
    H1y2 = Si_up2[:,1][i]+H21_shift_from_Si[1]
    if H1y2 < 0:
        H1y2 += max(Si_y)+y_pbc_dist
    H1z2 = Si_up2[:,2][i]+H21_shift_from_Si[2]
    H1id2 = int(len(Atoms))+i+1+5*len(Si_up2)
    H21.append([H1id2, 7, 0.06, H1x2, H1y2, H1z2])
    
    H2x2 = Si_up2[:,0][i]+H22_shift_from_Si[0]
    H2y2 = Si_up2[:,1][i]+H22_shift_from_Si[1]
    if H2y2 < 0:
        H2y2 += max(Si_y)+y_pbc_dist
    H2z2 = Si_up2[:,2][i]+H22_shift_from_Si[2]
    H2id2 = int(len(Atoms))+i+1+6*len(Si_up2)
    H22.append([H2id2, 7, 0.06, H2x2, H2y2, H2z2])
    
    H3x2 = Si_up2[:,0][i]+H23_shift_from_Si[0]
    H3y2 = Si_up2[:,1][i]+H23_shift_from_Si[1]
    if H3y2 < 0:
        H3y2 += max(Si_y)+y_pbc_dist
    H3z2 = Si_up2[:,2][i]+H23_shift_from_Si[2]
    H3id2 = int(len(Atoms))+i+1+7*len(Si_up2)
    H23.append([H3id2, 7, 0.06, H3x2, H3y2, H3z2])
    
    Bonds.append([1000, 2, int(Si_up2[:,3][i]), C2id])    # Si - C Numbering will do later
    Bonds.append([1000, 3, C2id, H1id2])    # C - H1
    Bonds.append([1000, 3, C2id, H2id2])    # C - H2
    Bonds.append([1000, 3, C2id, H3id2])    # C - H3
    
    Angles.append([1000, 2, Si_up2[:,3][i], C2id, H1id2])    #Si - C - H1
    Angles.append([1000, 2, Si_up2[:,3][i], C2id, H2id2])    #Si - C - H2
    Angles.append([1000, 2, Si_up2[:,3][i], C2id, H3id2])    #Si - C - H3
    Angles.append([1000, 3, H1id2, C2id, H2id2])    #H1 - C - H2
    Angles.append([1000, 3, H1id2, C2id, H3id2])    #H1 - C - H3
    Angles.append([1000, 3, H2id2, C2id, H3id2])    #H2 - C - H3
    Angles.append([1000, 4, Cid, Si_up2[:,3][i], C2id])    #C - Si - C
    
    
    


Atoms = Atoms + C1 + H11 + H12 + H13 + C2 + H21 + H22 + H23
print('Add CH3: '+str(sum(array(Atoms)[:,2])))

#--------- dublicate the whole slub with rotation around Y-axis------------
max_At_temp = max(array(Atoms)[:,5])


Bonds1 = array(Bonds)
for i in range(len(Bonds)):
    Bonds.append([1000, Bonds1[:,1][i], Bonds1[:,2][i]+len(Atoms), Bonds1[:,3][i]+len(Atoms)])
Angles1 = array(Angles)
for i in range(len(Angles)):
    Angles.append([1000, Angles1[:,1][i], Angles1[:,2][i]+len(Atoms), Angles1[:,3][i]+len(Atoms), Angles1[:,4][i]+len(Atoms)])




# Define the rotation matrix for a Y-axis rotation
angle_radians = np.radians(180)
rotation_matrix = np.array([[np.cos(angle_radians), 0, np.sin(angle_radians)], [0, 1, 0], [-np.sin(angle_radians), 0, np.cos(angle_radians)]])
Atoms_len = len(Atoms)
Atoms1 = array(Atoms)
for i in range(len(Atoms)):
    new_coord = np.dot(rotation_matrix, array([Atoms1[:,3][i], Atoms1[:,4][i], Atoms1[:,5][i]]))
    X = new_coord[0]+60.33338400000002
    Y = new_coord[1]
    Z = new_coord[2]+58.9266 #94.3022
    
    Atoms.append([Atoms1[:,0][i]+Atoms_len, Atoms1[:,1][i], Atoms1[:,2][i], X, Y, Z])
    

print('dublicated: '+str(sum(array(Atoms)[:,2])))


#--------------------------------------------------- WATER ------------------------

first_atom_row_H2O = 36

data_atoms_H2O = loadtxt(water_input, skiprows = row_num_H2O('Atoms #full')+2, max_rows = row_num_H2O('Bonds #full')-row_num_H2O('Atoms #full')-3)
H2O_id = data_atoms_H2O[:,0]
H2O_mol = data_atoms_H2O[:,1]
H2O_type = data_atoms_H2O[:,2]
H2O_charge = data_atoms_H2O[:,3]
H2O_x = data_atoms_H2O[:,4]
H2O_y = data_atoms_H2O[:,5]
H2O_z = data_atoms_H2O[:,6]

data_bonds_H2O = loadtxt(water_input, skiprows = row_num_H2O('Bonds #full')+2, max_rows = row_num_H2O('Angles #full')-row_num_H2O('Bonds #full')-2)
B_ind1_H2O = data_bonds_H2O[:,2]
B_ind2_H2O = data_bonds_H2O[:,3]

data_angs_H2O = loadtxt(water_input, skiprows = row_num_H2O('Angles #full')+2)
ANG_ind1_H2O = data_angs_H2O[:,2]
ANG_ind2_H2O = data_angs_H2O[:,3]
ANG_ind3_H2O = data_angs_H2O[:,4]


#move H2O
H2O_x_shift = 0
H2O_y_shift = -min(H2O_y)+((max(Si_y)+y_pbc_dist)-(max(H2O_y)-min(H2O_y)))/2
H2O_z_shift = max_At_temp+2-min(H2O_z)


len_at = len(Atoms)
for i in range(len(H2O_x)):
    Atoms.append([H2O_id[i]+len_at, int(H2O_type[i]+7), H2O_charge[i], H2O_x[i]+H2O_x_shift, H2O_y[i]+H2O_y_shift, H2O_z[i]+H2O_z_shift])
for i in range(len(B_ind1_H2O)):
    Bonds.append([1000, 4, B_ind1_H2O[i]+len_at, B_ind2_H2O[i]+len_at])
for i in range(len(ANG_ind1_H2O)):
    Angles.append([1000, 5, ANG_ind1_H2O[i]+len_at, ANG_ind2_H2O[i]+len_at, ANG_ind3_H2O[i]+len_at])

#-----------------------------------------------------------------------------------

print('Add water: '+str(sum(array(Atoms)[:,2])))

Bonds = array(Bonds)
Angles = array(Angles)
Atoms = array(Atoms)




for i in range(len(Bonds)):
    Bonds[:,0][i] = int(i+1)
for i in range(len(Angles)):
    Angles[:,0][i] = int(i+1)
for i in range(len(Atoms)):
    if Atoms[:,0][i] in Si_up2[:,3]:
        Atoms[:,1][i] = int(5)
        Atoms[:,2][i] = 1.05
    Atoms[:,5][i] +=  -O1down_shift_from_Si[2]
        





#-----------------------------WRITING------------------------------------

with open(output_file, 'w') as fdata:
    fdata.write('LAMMPS data file. CGCMM style. atom_style full\n')
    
    fdata.write('{} atoms\n'.format(str(len(Atoms))))
    fdata.write('{} bonds\n'.format(str(len(Bonds))))
    fdata.write('{} angles\n'.format(len(Angles)))   
    fdata.write('\n9 atom types\n4 bond types\n5 angle types\n\n')
   
    #System dymentions section
    
    fdata.write('{} {} xlo xhi\n'.format(min(Atoms[:,3]), max(Si_x)+x_pbc_dist))
    fdata.write('{} {} ylo yhi\n'.format(min(Atoms[:,4]), max(Si_y)+y_pbc_dist)) 
    fdata.write('{} {} zlo zhi\n'.format(-50, max(array(Atoms)[:,5])+50)) 
    
    fdata.write('\n')
    
    fdata.write('Masses\n\n')
    
    fdata.write('1 28.0855 # Si\n')
    fdata.write('2 15.999000 # O \n')
    fdata.write('3 15.999000 # O \n')
    fdata.write('4 1.00784   # H\n')
    fdata.write('5 28.0855 # Si\n')
    fdata.write('6 12.011 # C\n')
    fdata.write('7 1.00784   # H\n')
    fdata.write('8 15.999000 # O \n')
    fdata.write('9 1.00784   # H\n\n')
    
    #Atoms section
    fdata.write('Atoms #full\n\n')

    #SiO2
    i =0
    for pos in range(len(Atoms)):
        fdata.write(str(int(Atoms[:,0][pos]))+' '+str(1)+' '+str(int(Atoms[:,1][pos]))+' '+str(Atoms[:,2][pos])+' '+str(Atoms[:,3][pos])+' '+str(Atoms[:,4][pos])+' '+str(Atoms[:,5][pos])+'\n')
        i += 1
    N_last_Si_atom = i+1
    i = 0
    



    #Bonds section  
    fdata.write('\nBonds #full\n\n')    

    #SiO2
    i=0
    for pos in range(len(Bonds)):
        fdata.write(str(int(Bonds[:,0][pos]))+' '+str(int(Bonds[:,1][pos]))+' '+str(int(Bonds[:,2][pos]))+' '+str(int(Bonds[:,3][pos]))+'\n')
        i += 1
    N_last_Si_bond = i+1
    

    
    #Angles section
    fdata.write('\nAngles #full\n\n')   
    
    #SiO2
    i=0
    for pos in range(len(Angles)):
        fdata.write(str(int(Angles[:,0][pos]))+' '+str(int(Angles[:,1][pos]))+' '+str(int(Angles[:,2][pos]))+' '+str(int(Angles[:,3][pos]))+' '+str(int(Angles[:,4][pos]))+'\n')
        i += 1
    N_last_Si_angle = i+1
    



fdata.close()



    