#################################################################################################
####Addapted script from the extract.py script to be able to look at a forced angle scan where I scan the bond dist and constrain the bond angle over a range of angles. Pulls copper and scanned oxygen distance as well as the energy, then converts energy from Hartrees to kcal/mol. Does not save any new files or write xyz coords out anywhere. Compares energy from scans to a Potential
###############################################################################################################


#PREAMBLE:	

import sys
import math
import numpy as np
import matplotlib.pyplot as plt



#QM_energy is the energy from qm computed for a given frame/dist  
QM_energy1 = []
QM_energy2 = []
QM_energy3 = []
QM_energy4 = []
QM_energy5 = []

#cu_o_bond_dist is the list of coordination distances between copper and the oxygen being incrementally moved closer to the copper 1 -> 18 corresponding to increasing angle (90-180)
cu_o_bond_dist1 = []
cu_o_bond_dist2 = []
cu_o_bond_dist3 = []
cu_o_bond_dist4 = []
cu_o_bond_dist5 = []

#sys.argv[0] is extract script and the sys.argv after is the different .log files for angles in increasing order (90 ->180). The name_split splits the name of the .log file and the angle_scanned is the angle the the .log file is named with.
log_file1 = sys.argv[1]
name_split1 = ' '.join(' '.join(log_file1.split('_')).split('.')).split()
angle_scanned1 = name_split1[4]

log_file2 = sys.argv[2]
name_split2 = ' '.join(' '.join(log_file2.split('_')).split('.')).split()
angle_scanned2 = name_split2[4]

log_file3 = sys.argv[3]
name_split3 = ' '.join(' '.join(log_file3.split('_')).split('.')).split()
angle_scanned3 = name_split3[4]

log_file4 = sys.argv[4]
name_split4 = ' '.join(' '.join(log_file4.split('_')).split('.')).split()
angle_scanned4 = name_split4[4]

log_file5 = sys.argv[5]
name_split5 = ' '.join(' '.join(log_file5.split('_')).split('.')).split()
angle_scanned5 = name_split5[4]



#conv_factor: changes charges in subroutine to charges used in AMBER force fields. Gromacs and others use 18.22615. 
conv_factor = 18.2223


#atom_name1 and atom_number1 taken from first log file. 
atom_name1 = []
atom_name2 = []
atom_name3 = []
atom_name4 = []
atom_name5 = []
#atom_number refers to the order of the atoms. Thus if copper is always the first atom, it has a atom_number of 1.
atom_number1 = []
atom_number2 = []
atom_number3 = []
atom_number4 = []
atom_number5 = []
#cartesian values
x1 = []
x2 = []
x3 = []
x4 = []
x5 = []
y1 = []
y2 = []
y3 = []
y4 = []
y5 = []
z1 = []
z2 = []
z3 = []
z4 = []
z5 = []
#---------------------------------------------------------------------------------------------------
#Subroutines

def atomNumber2ElementName(atom_number):

        if atom_number == 1:
                element_name = "H"
        elif atom_number == 2:
                element_name = "He"
        elif atom_number == 6:
                element_name = "C"
        elif atom_number == 7:
                element_name = "N"
        elif atom_number == 8:
                element_name = "O"
        elif atom_number == 29:
                element_name = "Cu"

        return element_name

def MorsePotentialOFAllOxygens(coords, names):
        #De and re were values "measured" from the qm scan
        #1h2o
        #cu_o_De = 17.4015979 
        #cu_o_re = 1.9693 
        #cu_o_a = 3.200655 
	

	n_atoms = coords.shape[0]
	De_array0 = np.array([0., 17.4015979, 0.])
	De_array1 = np.array([17.4015979, 0., 0.])
	De_array2 = np.array([0., 0.,0.])
	De_array = np.concatenate(([De_array0], [De_array1], [De_array2]), axis=0)
	re_array0 = np.array([0., 1.9693, 0.])
	re_array1 = np.array([1.9693, 0., 0.])
	re_array2 = np.array([0., 0.,0.])
	re_array = np.concatenate(([re_array0], [re_array1], [re_array2]), axis=0)
	a_array0 = np.array([0., 3.200655, 0.])
	a_array1 = np.array([3.200655, 0., 0.])
	a_array2 = np.array([0., 0.,0.])
	a_array = np.concatenate(([a_array0], [a_array1], [a_array2]), axis=0)

        atom_type = np.zeros(n_atoms, dtype=int)
	for atom in range(n_atoms):
		if str(names[atom]) == "Cu":
			atom_type[atom] = 0
		elif str(names[atom]) == "O":
			atom_type[atom] = 1
		elif str(names[atom]) == "H":
			atom_type[atom] = 2

	potential = 0
        for atom1 in range(n_atoms - 1):
                for atom2 in range(atom1+1, n_atoms):
                        dist_vec = coords[atom1,:] - coords[atom2,:]
                        dist2 = np.dot(dist_vec,dist_vec)
                        dist = math.sqrt(dist2)
                     	typ1 = atom_type[atom1]
			typ2 = atom_type[atom2]
			De = De_array[typ1, typ2]
			re = re_array[typ1, typ2]
			a = a_array[typ1, typ2]
                        Morse = (De * (1 - math.exp(-a*(dist - re)))**2)
                        potential += float(Morse)

        return potential









#------------------------------------------------------------------------------------------
#MAIN PROGRAM:


with open(log_file1) as input:
	switch = 0
	count = 0 
	count2 = 0	
	#initial_parameters is used to identify where the "Scan" term is being found to help identify which oxygen is the one being used for the scan
	initial_parameters = 0
	for line in input: 
		d=line.split()
		if len(d)>1:
			#finding the "Scan" term in the log file
			if d[0] == "!" and d[1] == "Initial" and d[2] == "Parameters":
				initial_parameters += 1
			#finding which R is assigned for the copper-oxygen bond distance
			if initial_parameters == 1:
				if d[0] == "!" and d[1] == "R1" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R2" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R3" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R4" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R5" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
			if d[0] == "NAtoms=":
				n_atoms = int(d[1])
			if d[0]=="Optimization" and d[1]=="completed.": 
				switch += 1
			if switch == 1:
				if important_oxygen_name == d[1] and important_oxygen_definition == d[2]:
					#print d[3]
					cu_o_bond_dist1.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy1.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
					if count2 > 3 and count2 < n_atoms + 4:
	                                        #print d[0] 
        	                                atom_number1.append(int(d[0]))
                	                        atom_name1.append(atomNumber2ElementName(int(d[1])))
						x1.append(float(d[3]))
                                                y1.append(float(d[4]))
                                                z1.append(float(d[5]))

				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0


with open(log_file2) as input:
	switch = 0
	count = 0 
	count2 = 0	
	#initial_parameters is used to identify where the "Scan" term is being found to help identify which oxygen is the one being used for the scan
	initial_parameters = 0
	for line in input: 
		d=line.split()
		if len(d)>1:
			#finding the "Scan" term in the log file
			if d[0] == "!" and d[1] == "Initial" and d[2] == "Parameters":
				initial_parameters += 1
			#finding which R is assigned for the copper-oxygen bond distance
			if initial_parameters == 1:
				if d[0] == "!" and d[1] == "R1" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R2" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R3" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R4" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R5" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
			if d[0] == "NAtoms=":
				n_atoms = int(d[1])
			if d[0]=="Optimization" and d[1]=="completed.": 
				switch += 1
			if switch == 1:
				if important_oxygen_name == d[1] and important_oxygen_definition == d[2]:
					#print d[3]
					cu_o_bond_dist2.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy2.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
                                        if count2 > 3 and count2 < n_atoms + 4:
                                                #print d[0] 
                                                atom_number2.append(int(d[0]))
                                                atom_name2.append(atomNumber2ElementName(int(d[1])))
                                                x2.append(float(d[3]))
                                                y2.append(float(d[4]))
                                                z2.append(float(d[5]))

				if count2 == n_atoms + 70:
					switch = 0
					count = 0
					count2 = 0 


with open(log_file3) as input:
	switch = 0
	count = 0 
	count2 = 0	
	#initial_parameters is used to identify where the "Scan" term is being found to help identify which oxygen is the one being used for the scan
	initial_parameters = 0
	for line in input: 
		d=line.split()
		if len(d)>1:
			#finding the "Scan" term in the log file
			if d[0] == "!" and d[1] == "Initial" and d[2] == "Parameters":
				initial_parameters += 1
			#finding which R is assigned for the copper-oxygen bond distance
			if initial_parameters == 1:
				if d[0] == "!" and d[1] == "R1" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R2" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R3" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R4" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R5" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
			if d[0] == "NAtoms=":
				n_atoms = int(d[1])
			if d[0]=="Optimization" and d[1]=="completed.": 
				switch += 1
			if switch == 1:
				if important_oxygen_name == d[1] and important_oxygen_definition == d[2]:
					#print d[3]
					cu_o_bond_dist3.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy3.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
                                        if count2 > 3 and count2 < n_atoms + 4:
                                                #print d[0] 
                                                atom_number3.append(int(d[0]))
                                                atom_name3.append(atomNumber2ElementName(int(d[1])))
                                                x3.append(float(d[3]))
                                                y3.append(float(d[4]))
                                                z3.append(float(d[5]))

				if count2 == n_atoms + 70:
					switch = 0
					count = 0 
					count2 = 0


with open(log_file4) as input:
	switch = 0
	count = 0 
	count2 = 0	
	#initial_parameters is used to identify where the "Scan" term is being found to help identify which oxygen is the one being used for the scan
	initial_parameters = 0
	for line in input: 
		d=line.split()
		if len(d)>1:
			#finding the "Scan" term in the log file
			if d[0] == "!" and d[1] == "Initial" and d[2] == "Parameters":
				initial_parameters += 1
			#finding which R is assigned for the copper-oxygen bond distance
			if initial_parameters == 1:
				if d[0] == "!" and d[1] == "R1" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R2" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R3" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R4" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R5" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
			if d[0] == "NAtoms=":
				n_atoms = int(d[1])
			if d[0]=="Optimization" and d[1]=="completed.": 
				switch += 1
			if switch == 1:
				if important_oxygen_name == d[1] and important_oxygen_definition == d[2]:
					#print d[3]
					cu_o_bond_dist4.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy4.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
                                        if count2 > 3 and count2 < n_atoms + 4:
                                                #print d[0] 
                                                atom_number4.append(int(d[0]))
                                                atom_name4.append(atomNumber2ElementName(int(d[1])))
                                                x4.append(float(d[3]))
                                                y4.append(float(d[4]))
                                                z4.append(float(d[5]))

				if count2 == n_atoms + 70:
					switch = 0
					count = 0 
					count2 = 0


with open(log_file5) as input:
	switch = 0
	count = 0 
	count2 = 0	
	#initial_parameters is used to identify where the "Scan" term is being found to help identify which oxygen is the one being used for the scan
	initial_parameters = 0
	for line in input: 
		d=line.split()
		if len(d)>1:
			#finding the "Scan" term in the log file
			if d[0] == "!" and d[1] == "Initial" and d[2] == "Parameters":
				initial_parameters += 1
			#finding which R is assigned for the copper-oxygen bond distance
			if initial_parameters == 1:
				if d[0] == "!" and d[1] == "R1" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R2" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R3" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R4" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
				elif d[0] == "!" and d[1] == "R5" and d[4] == "Scan":
					important_oxygen_name = d[1]
					important_oxygen_definition = d[2]
			if d[0] == "NAtoms=":
				n_atoms = int(d[1])
			if d[0]=="Optimization" and d[1]=="completed.": 
				switch += 1
			if switch == 1:
				if important_oxygen_name == d[1] and important_oxygen_definition == d[2]:
					#print d[3]
					cu_o_bond_dist5.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy5.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
                                        if count2 > 3 and count2 < n_atoms + 4:
                                                #print d[0] 
                                                atom_number5.append(int(d[0]))
                                                atom_name5.append(atomNumber2ElementName(int(d[1])))
                                                x5.append(float(d[3]))
                                                y5.append(float(d[4]))
                                                z5.append(float(d[5]))

				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0


#------Fitted LJ------
#LJ(9-6) A,B,C values from 1 water system
A = 5038.05204605 
B = 1326.20662767 
C = -1076790.05551
#sigma = (A/B)**(-3)
#epsilon = (B**2) / (4*A)

Fit_LJ = []
for r in cu_o_bond_dist1[1:26]:
        l = A/(r**9) + B/(r**6) - C  
        Fit_LJ.append(l)
#element1inFitLJ = Fit_LJ[0]
#for i in range(len(Fit_LJ)):
#        Fit_LJ[i] -= element1inFitLJ






#zero out the relative qm energies
first_element_in_qm_energy1 = QM_energy1[0]
for i in range(len(QM_energy1)):
	QM_energy1[i] -= first_element_in_qm_energy1

first_element_in_qm_energy2 = QM_energy2[0]
for i in range(len(QM_energy2)):
	QM_energy2[i] -= first_element_in_qm_energy2

first_element_in_qm_energy3 = QM_energy3[0]
for i in range(len(QM_energy3)):
	QM_energy3[i] -= first_element_in_qm_energy3

first_element_in_qm_energy4 = QM_energy4[0]
for i in range(len(QM_energy4)):
	QM_energy4[i] -= first_element_in_qm_energy4

first_element_in_qm_energy5 = QM_energy5[0]
for i in range(len(QM_energy5)):
	QM_energy5[i] -= first_element_in_qm_energy5

n_frames = (len(x1) / n_atoms) -1

x1 = np.asarray(x1)
y1 = np.asarray(y1)
z1 = np.asarray(z1)
x2 = np.asarray(x2)
y2 = np.asarray(y2)
z2 = np.asarray(z2)
x3 = np.asarray(x3)
y3 = np.asarray(y3)
z3 = np.asarray(z3)
x4 = np.asarray(x4)
y4 = np.asarray(y4)
z4 = np.asarray(z4)
x5 = np.asarray(x5)
y5 = np.asarray(y5)
z5 = np.asarray(z5)


#declare coordinate matrix to contain all positions for every step
xyz1 = np.empty((n_frames,n_atoms,3),dtype=float)
xyz2 = np.empty((n_frames,n_atoms,3),dtype=float)
xyz3 = np.empty((n_frames,n_atoms,3),dtype=float)
xyz4 = np.empty((n_frames,n_atoms,3),dtype=float)
xyz5 = np.empty((n_frames,n_atoms,3),dtype=float)

# populate coordinate matrix
line_count1 = 0
for frame in range(n_frames):
        for atom in range(n_atoms):
                xyz1[frame,atom,0] = x1[line_count1]
                xyz1[frame,atom,1] = y1[line_count1]
                xyz1[frame,atom,2] = z1[line_count1]
                line_count1 += 1

# populate coordinate matrix
line_count2 = 0
for frame in range(n_frames):
        for atom in range(n_atoms):
                xyz2[frame,atom,0] = x2[line_count2]
                xyz2[frame,atom,1] = y2[line_count2]
                xyz2[frame,atom,2] = z2[line_count2]
                line_count2 += 1

# populate coordinate matrix
line_count3 = 0
for frame in range(n_frames):
        for atom in range(n_atoms):
                xyz3[frame,atom,0] = x3[line_count3]
                xyz3[frame,atom,1] = y3[line_count3]
                xyz3[frame,atom,2] = z3[line_count3]
                line_count3 += 1

# populate coordinate matrix
line_count4 = 0
for frame in range(n_frames):
        for atom in range(n_atoms):
                xyz4[frame,atom,0] = x4[line_count4]
                xyz4[frame,atom,1] = y4[line_count4]
                xyz4[frame,atom,2] = z4[line_count4]
                line_count4 += 1

# populate coordinate matrix
line_count5 = 0
for frame in range(n_frames):
        for atom in range(n_atoms):
                xyz5[frame,atom,0] = x5[line_count5]
                xyz5[frame,atom,1] = y5[line_count5]
                xyz5[frame,atom,2] = z5[line_count5]
                line_count5 += 1

M_all_oxygens1 = []
M_all_oxygens2 = []
M_all_oxygens3 = []
M_all_oxygens4 = []
M_all_oxygens5 = []

for frame in range(n_frames):
        potential1 = MorsePotentialOFAllOxygens(xyz1[frame,:,:], atom_name1[0:n_atoms])
        M_all_oxygens1.append(potential1)
#zero out the morse potential
first_element_in_Morse1 = M_all_oxygens1[0]
for i in range(len(M_all_oxygens1)):
        M_all_oxygens1[i] -= first_element_in_Morse1

for frame in range(n_frames):
        potential2 = MorsePotentialOFAllOxygens(xyz2[frame,:,:], atom_name2[0:n_atoms])
        M_all_oxygens2.append(potential2)
first_element_in_Morse2 = M_all_oxygens2[0]
for i in range(len(M_all_oxygens2)):
        M_all_oxygens2[i] -= first_element_in_Morse2

for frame in range(n_frames):
        potential3 = MorsePotentialOFAllOxygens(xyz3[frame,:,:], atom_name3[0:n_atoms])
        M_all_oxygens3.append(potential3)
first_element_in_Morse3 = M_all_oxygens3[0]
for i in range(len(M_all_oxygens3)):
        M_all_oxygens3[i] -= first_element_in_Morse3

for frame in range(n_frames):
        potential4 = MorsePotentialOFAllOxygens(xyz4[frame,:,:], atom_name4[0:n_atoms])
        M_all_oxygens4.append(potential4)
first_element_in_Morse4 = M_all_oxygens4[0]
for i in range(len(M_all_oxygens4)):
        M_all_oxygens4[i] -= first_element_in_Morse4

for frame in range(n_frames):
        potential5 = MorsePotentialOFAllOxygens(xyz5[frame,:,:], atom_name5[0:n_atoms])
        M_all_oxygens5.append(potential5)
first_element_in_Morse5 = M_all_oxygens5[0]
for i in range(len(M_all_oxygens5)):
        M_all_oxygens5[i] -= first_element_in_Morse5




#-------------------------------------------------------------------------------
#OTHER UTILITIES/PLOTTING:

plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.plot(cu_o_bond_dist1[1:26], QM_energy1[0:25], "g^")
#plt.plot(cu_o_bond_dist2[1:26], QM_energy2[0:25], "b^")
#plt.plot(cu_o_bond_dist3[1:26], QM_energy3[0:25], "r^")
#plt.plot(cu_o_bond_dist4[1:26], QM_energy4[0:25], "k^")
#plt.plot(cu_o_bond_dist5[1:26], QM_energy5[0:25], "m^")
plt.plot(cu_o_bond_dist1[0:23], M_all_oxygens1[0:23], "gs")
#plt.plot(cu_o_bond_dist2[0:24], M_all_oxygens2[0:24], "bs")
#plt.plot(cu_o_bond_dist3[0:24], M_all_oxygens3[0:24], "rs")
#plt.plot(cu_o_bond_dist4[0:24], M_all_oxygens4[0:24], "ks")
#plt.plot(cu_o_bond_dist5[0:24], M_all_oxygens5[0:24], "ms")
plt.plot(cu_o_bond_dist1[0:23], Fit_LJ[0:23],'g+' )


#plt.legend([angle_scanned5,'Morse of 180'], fontsize='10', loc='upper right', shadow =True, ncol=1, numpoints = 1)
plt.legend([angle_scanned1,'Morse of 100', 'Fitted LJ(9-6)'], fontsize='10', loc='upper right', shadow =True, ncol=1, numpoints = 1)
#plt.legend(['Fit LJ9-6', angle_scanned1, angle_scanned2,  angle_scanned3, angle_scanned4, angle_scanned5,'Morse of 100', 'Morse of 110', 'Morse of 150', 'Morse of 155', 'Morse of 180'], fontsize='10', loc='upper right', shadow =True, ncol=1, numpoints = 1)
#plt.legend([angle_scanned1, angle_scanned2,  angle_scanned3, angle_scanned4, angle_scanned5,'Morse of 100', 'Morse of 110', 'Morse of 150', 'Morse of 155', 'Morse of 180'], fontsize='10', loc='upper right', shadow =True, ncol=1, numpoints = 1)


plt.ylabel('Energy (kcal/mol)')
plt.xlabel('Copper-Oxygen Distance ($\AA$)')
plt.xlim((0,4))
plt.axhline(y = 0, color = "k")
plt.title(r'Angle Scan 100 2-Water System', size='14')
plt.suptitle(r'Morse Parameters from 1 Water System Bond Scan', size='14')
plt.show()







