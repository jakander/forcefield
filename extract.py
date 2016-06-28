
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
	#o_charge = -0.84 * conv_factor
	#h_charge = 0.42 * conv_factor
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


def MorsePotentialOFAllOxygens(coords, names):
	#De and re were values "measured" from the qm scan
	#2h2o
	#cu_o_De = 40.9896 / 2
	#cu_o_re = 1.9975 / 2 

	#6h2o
	#cu_o_De = 13.68250 / 2
	#cu_o_re = 2.3 / 2 

	#1h2o
	cu_o_De = 17.4015979 / 2
	cu_o_re = 1.9693 / 2 
	cu_o_a = 3.200655 / 2
	#cu_o_a = 2.99999999 / 2

	#a parameter found from AMBER94
	#cu_o_a = 2.6368 / 2

	#sigma is in units of angstroms
	#epsilon in units of kcal/mol
	o_sigma = 3.090
	o_epsilon = 0.0088629
	
	dist2_cut = 1.2 * 1.2
	cu_charge = 0.0 * conv_factor
	#Oxygen and Hydrogen charges are the charges used in the Tip3P solvent model
	#o_charge = -0.84 * conv_factor
	#h_charge = 0.42 * conv_factor
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
#			sigma = (atom_sigma[atom1] + atom_sigma[atom2]) / 2.0 
#			#geometric mean for epsilon 
#			epsilon = (atom_epsilon[atom1] * atom_epsilon[atom2]) ** 0.5
#			LJ = 4.0 * epsilon * ((sigma / dist) ** 12.0 - (sigma / dist) ** 6.0)
#			potential += float(LJ)
#	                if dist2 > dist2_cut:
#	       	 	       dist = math.sqrt(dist2)
#	       	  	       coul_energy = atom_charge[atom1] * atom_charge[atom2] / dist
#			       potential += float(coul_energy)	

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

#QM_energy is the energy quantum mechanically computed for a given frame 
QM_energy = []

#cu_o_bond_dist is the coordination distance between copper and the oxygen being incrementally moved closer to the copper
cu_o_bond_dist = []

coords = {}
#sys.argv[1] is the script extract.py
log_file = sys.argv[1]
name = ' '.join(' '.join(log_file.split("_")).split('.')).split() 
print name
print log_file
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

MorseScannedWaterOnly = []
for i in range(len(cu_o_bond_dist)):
	#1water 
	De = 17.4015979
	re = 1.9693 
	
	#6 water De and re
	#De = 13.68250
	#re = 2.3

	#2 water De and re
	#De = 40.340
	#re = 1.9975
	a = 3.200655
	Morse = (De * (1 - math.exp(-a*(cu_o_bond_dist[i] - re)))**2 - De)
	MorseScannedWaterOnly.append(float(Morse))


#compute the sigma and epsilon values for the Lennard Jones through a least squares approach
#repulsive_term = [7,8,9,10,11,12,13,14,15,16,17,18,19,20]
#colors_repulsive_term = ['bo','ro','ko','go','mo','bs','rs','ks','gs','ms','b^','k^','g^','m^']
repulsive_term = [9, 12]
colors_repulsive_term = ['ms','b^']
#repulsive_term = [15,16,17,18]
#colors_repulsive_term = ['gs','ms','b^','k^']
#repulsive_term = [13,14,15,16]
#colors_repulsive_term = ['rs','ks','gs','ms']
#repulsive_term = [10,11,12,13,14]
#colors_repulsive_term = ['go','mo','bs','rs','ks']
#repulsive_term = [9,10,11,12]
#colors_repulsive_term = ['ko','go','mo','bs']
#repulsive_term = [7,8,9,10]
#colors_repulsive_term = ['bo','ro','ko','go']
for m in range(len(repulsive_term)):
	r1 = []
	r2 = []
	r3 = []
	for i in cu_o_bond_dist[1:26]:
		l = i ** (-repulsive_term[m])
		k = -i ** (-6)
		r1.append(l)
		r2.append(k)
		r3.append(1)
	r_all = [r1, r2, r3]
	r_all = np.asmatrix(r_all)
	QM_energy = np.asarray(QM_energy)
	A, B , C= np.linalg.lstsq(r_all.T, QM_energy)[0]
	#residuals = np.linalg.lstsq(r_all.T,QM_energy)[1]
	#LJ(9-6) A,B,C values from 1 water system
	#A = 5038.05204605 
	#B = 1326.20662767 
	#C = -1076790.05551
	#LJ(20-6) A,B,C values from 1 water system
	#A = 57801.8233555 
	#B = 53.9734373032 
	#C = -1076793.99734
	#print A, B, C
	#print (m+7), "= ", residuals
	#sigma = (A/B)**(-3)
	#epsilon = (B**2) / (4*A)
	Fit_LJ = []
	for r in cu_o_bond_dist[1:26]:
		l = A/(r**repulsive_term[m]) + B/(r**6) - C
		Fit_LJ.append(l)
	element1inFitLJ = Fit_LJ[0]
	for i in range(len(Fit_LJ)):
		Fit_LJ[i] -= element1inFitLJ
	plt.plot(cu_o_bond_dist[1:26], Fit_LJ[0:25], colors_repulsive_term[m])
	
	
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


energy_list = []
for frame in range(n_frames):
	energy = computePotentialEnergy(xyz[frame,:,:],atom_name[0:n_atoms])
	energy_list.append(energy)

M_all_oxygens = []
for frame in range(n_frames):
	potential = MorsePotentialOFAllOxygens(xyz[frame,:,:], atom_name[0:n_atoms])
	M_all_oxygens.append(potential)



M_LJ_C = []
for i in energy_list:
	for j in MorseScannedWaterOnly:
		k = i + j
		M_LJ_C.append(k)




#zero out the morsescannedwater + LJ+C
first_element_in_M_LJ_C = M_LJ_C[0]
for i in range(len(M_LJ_C)):
	M_LJ_C[i] -= first_element_in_M_LJ_C

#zero out the morse potential
first_element_in_Morse = M_all_oxygens[0]
for i in range(len(M_all_oxygens)):
	M_all_oxygens[i] -= first_element_in_Morse

#zero out the relative potential energy
first_element_in_ff_pe = energy_list[0]
for i in range(len(energy_list)):
	energy_list[i] -= first_element_in_ff_pe

#zero out the relative qm energy
first_element_in_qm_energy = QM_energy[0]
for i in range(len(QM_energy)):
	QM_energy[i] -= first_element_in_qm_energy





plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
#plt.plot(cu_o_bond_dist[1:26], energy_list[0:25], "bo")
plt.plot(cu_o_bond_dist[1:25], QM_energy[0:24], "r^")
#plt.plot(cu_o_bond_dist[26:50], QM_energy[25:49], "g^")
#plt.plot(cu_o_bond_dist[51:75], QM_energy[50:74], "b^")
#plt.plot(cu_o_bond_dist[76:100], QM_energy[75:99], "k^")
#plt.plot(cu_o_bond_dist[101:125], QM_energy[100:124], "rs")
#plt.plot(cu_o_bond_dist[126:148], QM_energy[125:147], "gs")
#plt.plot(cu_o_bond_dist[1:26], Fit_LJ[0:25], "bo")
plt.plot(cu_o_bond_dist[1:26], M_all_oxygens[0:25], "gs")
#plt.plot(cu_o_bond_dist[1:26], M_PE[0:25], "ko")
#plt.plot(cu_o_bond_dist[1:26], M_LJ_C[0:25], "ko")
plt.ylabel('Energy (kcal/mol)')
plt.xlabel('Copper-Oxygen Distance ($\AA$)')
#plt.legend(['LJ + C', 'QM', 'Morse', 'M_PE'], fontsize='10', bbox_to_anchor=(.85, 0.865, 0.15, 0.0), loc=3, ncol=1, mode='expand', borderaxespad=0., numpoints = 1)
#plt.legend(['LJ15-6','LJ16-6','LJ17-6','LJ18-6','QM'], fontsize='10', bbox_to_anchor=(.83,.78, 0.17, 0.0), loc=3, ncol=1, mode='expand', borderaxespad=0., numpoints = 1)
#plt.legend(['LJ7-6','LJ8-6','LJ9-6','LJ10-6','LJ11-6','LJ12-6','LJ13-6','LJ14-6','LJ15-6','LJ16-6','LJ17-6','LJ18-6','LJ19-6','LJ20-6','QM'], fontsize='10', bbox_to_anchor=(.83, .373, 0.17, 0.0), loc=3, ncol=1, mode='expand', borderaxespad=0., numpoints = 1)
plt.legend(['LJ9-6', 'LJ12-6','QM','Morse'], fontsize='10', bbox_to_anchor=(.83, .78, 0.17, 0.0), loc=3, ncol=1, mode='expand', borderaxespad=0., numpoints = 1)
#plt.legend(['QM','LJ(16-6)'], fontsize='10', bbox_to_anchor=(.75, .9, 0.25, 0.0), loc=3, ncol=1, mode='expand', borderaxespad=0., numpoints = 1)
#plt.legend(['QM','LJ(16-6)'], fontsize='10', bbox_to_anchor=(.75, .9, 0.25, 0.0), loc=3, ncol=1, mode='expand', borderaxespad=0., numpoints = 1)
plt.xlim((0,4))
plt.axhline(y = 0, color = "k")
#plt.title(r'2D Coordination and Angle Potential Energy Scan', size='14')
#plt.suptitle(r'Copper and 2 Water', size='14')
#plt.suptitle(r'With Parameters from 1 Water System', size='14')
plt.title(r'Potential Energy Scan 1-Water System', size='14')
plt.show()

#print the QM energy and the scanned bond distances to new file to be read to determine the sigma and epsilon valuesfor Lennard Jones potentials 
#nf_paramLJ = open("paramLJ.xyz", "w")
#for i in range(len(QM_energy)):
#	nf_paramLJ.write("%4f\n" % (QM_energy))
#	nf_paramLJ.write("%4f\n" % (cu_o_bond_dist))
	
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


