
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

#conv_factor: changes charges in subroutine to charges used in AMBER force fields. Gromacs and others use 18.22615. 
conv_factor = 18.2223 

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


def computePotentialEnergy(coords, names):

	dist2_cut = 1.2 * 1.2
	cu_charge = 2.0 * conv_factor
	#Oxygen and Hydrogen charges are the charges used in the Tip3P solvent model
	#-.72 and .36 fit the qm scan a little better
	o_charge = -0.72 * conv_factor
	h_charge = 0.36 * conv_factor

	#sigma is in units of angstroms
	#epsilon in units of kcal/mol
	cu_sigma = 2.066 / (2.) ** (1. / 6.)
	o_sigma = 3.090
	cu_epsilon = 0.0427
	o_epsilon = 0.088629

	
	n_atoms = coords.shape[0]
	atom_sigma = np.zeros(n_atoms, dtype=float)
	for atom in range(n_atoms):
		if str(names[atom]) == "H":
			atom_sigma[atom] = 0.0
		elif str(names[atom]) == "O":
			atom_sigma[atom] = o_sigma
		elif str(names[atom]) == "Cu":
			atom_sigma[atom] = cu_sigma

	atom_epsilon = np.zeros(n_atoms, dtype=float)
	for atom in range(n_atoms):
		if names[atom] == "H":
			atom_epsilon[atom] = 0.0
		elif names[atom] == "O":
			atom_epsilon[atom] = o_epsilon
		elif names[atom] == "Cu":
			atom_epsilon[atom] = cu_epsilon

	atom_charge = np.empty(n_atoms,dtype=float)
	for atom in range(n_atoms):
		if names[atom] == "H":
			atom_charge[atom] = h_charge
		elif names[atom] == "O":
			atom_charge[atom] = o_charge
		elif names[atom] == "Cu":
			atom_charge[atom] = cu_charge


	energy = 0
	potential = 0
	for atom1 in range(n_atoms-1):

		for atom2 in range(atom1+1, n_atoms):
			dist_vec = coords[atom1,:]-coords[atom2,:]
			dist2 = np.dot(dist_vec,dist_vec)
			dist = math.sqrt(dist2)
			#arithemtic mean for sigma
			sigma = (atom_sigma[atom1] + atom_sigma[atom2]) / 2.0 
			#geometric mean for epsilon 
			epsilon = (atom_epsilon[atom1] * atom_epsilon[atom2]) ** 0.5
			LJ = 4.0 * epsilon * ((sigma / dist) ** 12.0 - (sigma / dist) ** 6.0)
			energy += float(LJ)
			if dist2 > dist2_cut:
				dist = math.sqrt(dist2)
				coul_energy = atom_charge[atom1] * atom_charge[atom2] / dist
				energy += float(coul_energy)	
	return energy 


def MorsePotential(coords, names):
	#De and re were values "measured" from the qm scan
	#2h2o
	cu_o_De = 40.9896 / 2
	#6h2o
	#cu_o_De = 14.1667 / 2
	cu_o_re = 1.9975 / 2 
	#a parameter found from AMBER94
	cu_o_a = 2.6368 / 2

	#sigma is in units of angstroms
	#epsilon in units of kcal/mol
	o_sigma = 3.090
	o_epsilon = 0.0088629
	
	dist2_cut = 1.2 * 1.2
	cu_charge = 0.0 * conv_factor
	#Oxygen and Hydrogen charges are the charges used in the Tip3P solvent model
	#-.72 and .36 fit the qm scan a little better
	o_charge = -0.72 * conv_factor
	h_charge = 0.0


	n_atoms = coords.shape[0]
	atom_De = np.zeros(n_atoms, dtype=float)
	for atom in range(n_atoms):
		if str(names[atom]) == "H":
			atom_De[atom] = 0.0
		elif str(names[atom]) == "O":
			atom_De[atom] = cu_o_De
		elif str(names[atom]) == "Cu":
			atom_De[atom] = cu_o_De
	
	atom_re = np.zeros(n_atoms, dtype=float)
	for atom in range(n_atoms):
		if str(names[atom]) == "H":
			atom_re[atom] = 0.0
		elif str(names[atom]) == "O":
			atom_re[atom] = cu_o_re
		elif str(names[atom]) == "Cu":
			atom_re[atom] = cu_o_re
	
	atom_a = np.zeros(n_atoms, dtype=float)
	for atom in range(n_atoms):
		if str(names[atom]) == "H":
			atom_a[atom] = 0.0
		elif str(names[atom]) == "O":
			atom_a[atom] = cu_o_a
		elif str(names[atom]) == "Cu":
			atom_a[atom] = cu_o_a
	
	atom_sigma = np.zeros(n_atoms, dtype=float)
	for atom in range(n_atoms):
		if str(names[atom]) == "H":
			atom_sigma[atom] = 0.0
		elif str(names[atom]) == "O":
			atom_sigma[atom] = o_sigma
		elif str(names[atom]) == "Cu":
			atom_sigma[atom] = 0.0

	atom_epsilon = np.zeros(n_atoms, dtype=float)
	for atom in range(n_atoms):
		if names[atom] == "H":
			atom_epsilon[atom] = 0.0
		elif names[atom] == "O":
			atom_epsilon[atom] = o_epsilon
		elif names[atom] == "Cu":
			atom_epsilon[atom] = 0.0
	
	atom_charge = np.empty(n_atoms,dtype=float)
	for atom in range(n_atoms):
		if names[atom] == "H":
			atom_charge[atom] = 0.0
		elif names[atom] == "O":
			atom_charge[atom] = o_charge
		elif names[atom] == "Cu":
			atom_charge[atom] = 0.0


	potential = 0
	
#	for i in range(len(cu_o_bond_dist[0:25])):
#		De = 40.9896
#		re = 1.9975
#		a = -2.6365
#		Morse = (De * (1 - math.exp(a*(cu_o_bond_dist[i] - re)))**2)
#		potential += float(Morse)
#	return potential
	for atom1 in range(n_atoms-1):
		for atom2 in range(atom1+1, n_atoms):
			dist_vec = coords[atom1,:] - coords[atom2,:]
			dist2 = np.dot(dist_vec,dist_vec)
			dist = math.sqrt(dist2)
			De = (atom_De[atom1] + atom_De[atom2])
			re = (atom_re[atom1] + atom_re[atom2])
			a = (atom_a[atom1] + atom_a[atom2])
			Morse = (De * (1 - math.exp(-a*(dist - re)))**2)
			potential += float(Morse)
			#arithemtic mean for sigma
			sigma = (atom_sigma[atom1] + atom_sigma[atom2]) / 2.0 
			#geometric mean for epsilon 
			epsilon = (atom_epsilon[atom1] * atom_epsilon[atom2]) ** 0.5
			LJ = 4.0 * epsilon * ((sigma / dist) ** 12.0 - (sigma / dist) ** 6.0)
			potential += float(LJ)
			if dist2 > dist2_cut:
				dist = math.sqrt(dist2)
				coul_energy = atom_charge[atom1] * atom_charge[atom2] / dist
				potential += float(coul_energy)	

	return potential



switch = 0
count = 0 
count2 = 0

#initial_parameters is used to identify where the "Scan" term is being found to help identify which oxygen is the one being used for the scan
initial_parameters = 0


#atom_number refers to the order of the atoms. Thus if copper is always the first atom, it has a atom_number of 1.
atom_number = []

#x, y, and z refer to the cartesian coordinate values for all of the atoms at each given frame
x = []
y = []
z = []
atom_name = []

#QM_energy is the quantum mechanically computed for a given frame 
QM_energy = []

#cu_o_bond_dist is the coordination distance between copper and the oxygen being incrementally moved closer to the copper
cu_o_bond_dist = []

coords = {}
#sys.argv[1] is the script extract.py
log_file = sys.argv[1]
with open(log_file) as input:
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
					cu_o_bond_dist.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
					if count2 > 3 and count2 < n_atoms + 4:
					#	print d[0] 
						atom_number.append(int(d[0]))
						atom_name.append(atomNumber2ElementName(int(d[1])))
						x.append(float(d[3]))
						y.append(float(d[4]))
						z.append(float(d[5]))
					if count2 == n_atoms + 70:
						switch = 0
						count = 0 
						count2 = 0

n_frames = len(x) / n_atoms
n_frames = n_frames - 1

x = np.asarray(x)
y = np.asarray(y)
z = np.asarray(z)

#declare coordinate matrix to contain all positions for every step
xyz = np.empty((n_frames,n_atoms,3),dtype=float)
# populate coordinate matrix
line_count = 0
for frame in range(n_frames):
	for atom in range(n_atoms):
		xyz[frame,atom,0] = x[line_count]
		xyz[frame,atom,1] = y[line_count]
		xyz[frame,atom,2] = z[line_count]
		line_count += 1

#print cu_o_bond_dist[0:25]


energy_list = []
for frame in range(n_frames):
	energy = computePotentialEnergy(xyz[frame,:,:],atom_name[0:n_atoms])
	energy_list.append(energy)

Morse_potential=[]
for frame in range(n_frames):
	potential = MorsePotential(xyz[frame,:,:], atom_name[0:n_atoms])
	Morse_potential.append(potential)

#zero out the morse potential
first_element_in_Morse = Morse_potential[0]
for i in range(len(Morse_potential)):
	Morse_potential[i] -= first_element_in_Morse

#zero out the relative potential energy
first_element_in_ff_pe = energy_list[0]
for i in range(len(energy_list)):
	energy_list[i] -= first_element_in_ff_pe

#zero out the relative qm energy
first_element_in_qm_energy = QM_energy[0]
for i in range(len(QM_energy)):
	QM_energy[i] -= first_element_in_qm_energy



plt.plot(cu_o_bond_dist[0:25], energy_list[0:25], "bo")
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.plot(cu_o_bond_dist[0:25], QM_energy[0:25], "r^")
plt.plot(cu_o_bond_dist[0:25], Morse_potential[0:25], "gs")
plt.ylabel('Energy (kcal/mol)')
plt.xlabel('Copper-Oxygen Distance ($\AA$)')
plt.legend(['LJ + C', 'QM', 'Morse'], fontsize='10', bbox_to_anchor=(.85, 0.865, 0.15, 0.0), loc=3, ncol=1, mode='expand', borderaxespad=0., numpoints = 1)
plt.xlim((0,4))
plt.axhline(y = 0, color = "k")
plt.title(r'Copper-6 Equatorial Water Coordination Potential Energy Scan', size='14')
plt.show()


# print coordinates to a new file 
nf_coords = open("coords.xyz", "w")
step_count = 0
for i in range(0,len(x)):
	if i%n_atoms == 0:
		step_count += 1
		nf_coords.write("%4d\n" % (n_atoms))
		nf_coords.write("step %4d\n" % (step_count))
	nf_coords.write("%4i %4s  %10.6f  %10.6f  %10.6f \n" % (atom_number[i], atom_name[i], x[i], y[i], z[i])) 
for l in range(0,n_frames):
	if l == 0:
		nf_coords.write("%+15.15f \n" % (QM_energy[l]))
nf_coords.close()


