
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

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

def computePotentialEnergy(coords, names, cu_o_atom_numbers, frame):
	
	thresh = 0.0001
	dist2_cut = 1.2 * 1.2
	cu_charge = 2.0 * conv_factor
	o_charge = -0.84 * conv_factor
	h_charge = 0.42 * conv_factor

	#sigma in angstrom
	#epsilon in kcal/mol
	cu_sigma = 2.066 / (2.) ** (1. / 6.)
	o_sigma = 3.090
	cu_epsilon = 0.0427
	o_epsilon = 0.0088629
	
	n_atoms = coords.shape[0]
	atom_sigma = np.zeros(n_atoms, dtype=float)
	for atom in range(n_atoms):
	#	print str(names[atom])
		if str(names[atom]) == "H":
			atom_sigma[atom] = 0.0
		elif str(names[atom]) == "O":
			atom_sigma[atom] = o_sigma
		elif str(names[atom]) == "Cu":
			atom_sigma[atom] = cu_sigma
	n_atoms = coords.shape[0]
	atom_epsilon = np.zeros(n_atoms, dtype=float)
	for atom in range(n_atoms):
		if names[atom] == "H":
			atom_epsilon[atom] = 0.0
		elif names[atom] == "O":
			atom_epsilon[atom] = o_epsilon
		elif names[atom] == "Cu":
			atom_epsilon[atom] = cu_epsilon

	n_atoms = coords.shape[0]
	atom_charge = np.empty(n_atoms,dtype=float)
	for atom in range(n_atoms):
		if names[atom] == "H":
			atom_charge[atom] = h_charge
		elif names[atom] == "O":
			atom_charge[atom] = o_charge
		elif names[atom] == "Cu":
			atom_charge[atom] = cu_charge

	energy = 0
	cu_o_dist = 0
	for atom1 in range(n_atoms-1):

		for atom2 in range(atom1+1, n_atoms):
			dist_vec = coords[atom1,:]-coords[atom2,:]
			dist2 = np.dot(dist_vec,dist_vec)
			dist = math.sqrt(dist2)
			if frame == 0:
				if  (names[atom1] == "Cu" and names[atom2] == "O") or (names[atom1] == "O" and names[atom2] == "Cu"):
					if abs(dist-0.1) < thresh :
						cu_o_atom_numbers[0] = atom1
						cu_o_atom_numbers[1] = atom2
						cu_o_dist = dist
						print cu_o_dist
			elif atom1 == cu_o_atom_numbers[0] and atom2 == cu_o_atom_numbers[1]:
				cu_o_dist = dist
			
			#arithemtic mean for sigma
			sigma = (atom_sigma[atom1] + atom_sigma[atom2]) / 2.0
			#geometric mean for epsilon 
			epsilon = (atom_epsilon[atom1] * atom_epsilon[atom2]) ** 0.5
			LJ = 4.0 * epsilon * ((sigma / dist) ** 12.0 - (sigma / dist) ** 6.0)
	#		print float(LJ), "LJ", sigma, "sigma", epsilon, "epsilon", dist, "dist"
			energy += float(LJ)
			if dist2 > dist2_cut:
				dist = math.sqrt(dist2)
				coul_energy = atom_charge[atom1] * atom_charge[atom2] / dist
				energy += float(coul_energy)
	
	return energy, cu_o_dist
		

switch=0
count=0
count2=0
testcount=0
#atom_number refers to the order of the atoms. Thus if copper is always the first atom, it has a atom_number of 1.
atom_number = []
x=[]
y=[]
z=[]
atom_name = []
QM_energy = []
#n_atoms = 19
coords = {}
log_file = sys.argv[1]
with open(log_file) as input:
	for line in input: 
		d=line.split()
		if len(d)>1:
			if d[0] == "NAtoms=":
				n_atoms = int(d[1])
			if d[0]=="Optimization" and d[1]=="completed.": 
				switch += 1
			if switch == 1:
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
				#	print d[4]
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

#print atom_name
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
	energy = computePotentialEnergy(xyz[frame,:,:], atom_name[0:n_atoms])
	energy_list.append(energy)

#zero out the relative potential energy
first_element_in_ff_pe = energy_list[0]
for i in range(len(energy_list)):
	energy_list[i] -= first_element_in_ff_pe

#zero out the relative qm energy
first_element_in_qm_energy = QM_energy[0]
for i in range(len(QM_energy)):
	QM_energy[i] -= first_element_in_qm_energy

plt.plot(energy_list, "bo")
plt.plot(QM_energy, "ro")
plt.ylabel('Kcal/mol')
plt.xlabel('Frames')
plt.show()


# print coordinates
nf_coords = open("coords.xyz", "w")
step_count = 0
for i in range(0,len(x)):
#for i in range(0,len(x)) and l in range(0,n_frames):
	if i%n_atoms == 0:
#	if i%n_atoms%l == 0:
		step_count += 1
		nf_coords.write("%4d\n" % (n_atoms))
		nf_coords.write("step %4d\n" % (step_count))
	#	nf_coords.write("step %4d %+20.20f \n" % (step_count, QM_energy[l]))
	nf_coords.write("%4i %4s  %10.6f  %10.6f  %10.6f \n" % (atom_number[i], atom_name[i], x[i], y[i], z[i])) 
#	nf_coords.write("%4s  %10.6f  %10.6f  %10.6f %+20.20f \n" % (atom_name[i], x[i], y[i], z[i], QM_energy[i])) 
for l in range(0,n_frames):
	if l == 0:
		nf_coords.write("%+15.15f \n" % (QM_energy[l]))
nf_coords.close()

#print "frames:", len(x)/n_atoms

