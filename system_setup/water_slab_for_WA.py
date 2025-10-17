# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 10:38:54 2024

@author: mandrolk1
"""



from pylab import *
import numpy as np

###--------- WATER CREATION ---------------###


#"lattice parameters" for water
a = 3.1        #[Å]

OH_bond = 0.964
shift_H = (OH_bond*2**(1/2))/2
shift_H = round(shift_H, 5)



system_size_x = int(floor(44.90066497049302/a))  # [elementary cells]   105.48342200000004 - wide   44.90066497049302 - narrow
system_size_y = int(floor(43.4/a))
system_size_z = int(floor(43.4/a))



positions_O = []

for i in range(system_size_x):
    for j in range(system_size_y):
        for k in range(system_size_z):
            base_positions = np.array([i,j,k])
            cart_positions = np.inner(a, base_positions)
            
            positions_O.append(np.round(cart_positions + 0.5*a, 5))


positions_H = []                
for pos in positions_O:
    positions_H.append(pos+[0,shift_H,shift_H])
    positions_H.append(pos+[0,-shift_H,shift_H])




#---------to move all water molecules
x_shift = 0
y_shift = 0
z_shift = 0

# moving O
for i in range(len(positions_O)):
    positions_O[i][1] += y_shift     
    positions_O[i][2] += z_shift

# moving H    
for i in range(len(positions_H)):
    positions_H[i][1] += y_shift
    positions_H[i][2] += z_shift 
 #-----------------------------------------------   
    


#---------------------------------------------------------------------

with open('C:/LEMTA/Data/SiO2/ClayFF/creation/system_configurations/WA/water.data', 'w') as fdata:
    fdata.write('LAMMPS data file. CGCMM style. atom_style full\n')

    fdata.write('{} atoms\n'.format(len(positions_H)+len(positions_O)))
    fdata.write('{} bonds\n'.format(len(positions_H)))
    fdata.write('{} angles\n'.format(len(positions_O)))  
    fdata.write('{} dihedrals\n'.format(0))
    fdata.write('0 impropers\n2 atom types\n1 bond types\n1 angle types\n0 dihedral types\n0 improper types\n\n')
    
    #System dymentions section
    
    fdata.write('{} {} xlo xhi\n'.format(0, 108.6))
    fdata.write('{} {} ylo yhi\n'.format(0, 50))
    fdata.write('{} {} zlo zhi\n'.format(-10, 70)) 
    fdata.write('\n')
    
    #Coeffs indication section
    
    fdata.write('#Pair Coeffs\n#\n# 1 O\n# 2 H\n')
    fdata.write('#Bond Coeffs\n#\n# 1 O-H\n\n')
    fdata.write('#Angle Coeffs\n#\n# 1 H-O-H\n\n')

    
    #Masses section    
    
    fdata.write('Masses\n\n')
    fdata.write('1 15.999000 # O\n')
    fdata.write('2 1.0078400 # H\n\n')
                
    
    #Atoms section
    
    fdata.write('Atoms #full\n\n')
   



    #write atom positions
    prev_atom_num = 0                 
    prev_molecule_num = 0             
    a=0
    molecule_number = 0
    for i in range(len(positions_O)):
        fdata.write('{} {} 1 -0.8476 {} {} {}\n'.format(a+1+prev_atom_num, molecule_number+1+prev_molecule_num, *positions_O[i] ))
        fdata.write('{} {} 2 0.4238 {} {} {}\n'.format(a+2+prev_atom_num, molecule_number+1+prev_molecule_num, *positions_H[i*2]))
        fdata.write('{} {} 2 0.4238 {} {} {}\n'.format(a+3+prev_atom_num, molecule_number+1+prev_molecule_num, *positions_H[i*2+1]))
        a=a+3
        molecule_number+=1
    
        
    fdata.write('\n')  
    
        
    #write bonds
    fdata.write('Bonds #full\n\n')   
    a=0
    b=1
    for i in range(len(positions_O)):
        fdata.write('{} 1 {} {} \n'.format(b+1, a+1+prev_atom_num, a+2+prev_atom_num ))
        fdata.write('{} 1 {} {} \n'.format(b+2, a+1+prev_atom_num, a+3+prev_atom_num ))
        a=a+3
        b=b+2
    
        
    #write angles
    fdata.write('Angles #full\n\n')   
    prev_angle_num = 0
    a=1
    for i in range(len(positions_O)):
        fdata.write('{} 1 {} {} {}\n'.format(i+1+prev_angle_num, a+1+prev_atom_num, a+prev_atom_num, a+2+prev_atom_num ))
        a=a+3
    prev_angle_num = prev_angle_num+i+1
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
