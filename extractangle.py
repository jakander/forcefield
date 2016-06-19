######Addapted script from the extract.py script to be able to look at a forced angle scan where I scan the bond dist and constrain the bond angle over a range of angles. Pulls copper and scanned oxygen distance as well as the energy, then converts energy from Hartrees to kcal/mol. Does not save any new files or write xyz coords out anywhere.######



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
QM_energy6 = []
QM_energy7 = []
QM_energy8 = []
QM_energy9 = []
QM_energy10 = []
QM_energy11 = []
QM_energy12 = []
QM_energy13 = []
QM_energy14 = []
QM_energy15 = []
QM_energy16 = []
QM_energy17 = []
QM_energy18 = []

#cu_o_bond_dist is the list of coordination distances between copper and the oxygen being incrementally moved closer to the copper 1 -> 18 corresponding to increasing angle (90-180)
cu_o_bond_dist1 = []
cu_o_bond_dist2 = []
cu_o_bond_dist3 = []
cu_o_bond_dist4 = []
cu_o_bond_dist5 = []
cu_o_bond_dist6 = []
cu_o_bond_dist7 = []
cu_o_bond_dist8 = []
cu_o_bond_dist9 = []
cu_o_bond_dist10 = []
cu_o_bond_dist11 = []
cu_o_bond_dist12 = []
cu_o_bond_dist13 = []
cu_o_bond_dist14 = []
cu_o_bond_dist15 = []
cu_o_bond_dist16 = []
cu_o_bond_dist17 = []
cu_o_bond_dist18 = []

#sys.argv[0] is extract script and the sys.argv after is the different .log files for angles in increasing order (90 ->180). The name_split splits the name of the .log file and the angle_scanned is the angle the the .log file is named with.
log_file1 = sys.argv[1]]
name_split1 = ' '.join(' '.join(log_file1.split('_')).split('.')).split()
angle_scanned1 = name_split1[5]

log_file2 = sys.argv[2]
name_split2 = ' '.join(' '.join(log_file2.split('_')).split('.')).split()
angle_scanned2 = name_split2[5]

log_file3 = sys.argv[3]
name_split3 = ' '.join(' '.join(log_file3.split('_')).split('.')).split()
angle_scanned3 = name_split3[5]

log_file4 = sys.argv[4]
name_split4 = ' '.join(' '.join(log_file4.split('_')).split('.')).split()
angle_scanned4 = name_split4[5]

log_file5 = sys.argv[5]
name_split5 = ' '.join(' '.join(log_file5.split('_')).split('.')).split()
angle_scanned5 = name_split5[5]

log_file6 = sys.argv[6]
name_split6 = ' '.join(' '.join(log_file6.split('_')).split('.')).split()
angle_scanned6 = name_split6[5]

log_file7 = sys.argv[7]
name_split7 = ' '.join(' '.join(log_file7.split('_')).split('.')).split()
angle_scanned7 = name_split7[5]

log_file8 = sys.argv[8]
name_split8 = ' '.join(' '.join(log_file8.split('_')).split('.')).split()
angle_scanned8 = name_split8[5]

log_file9 = sys.argv[9]
name_split9 = ' '.join(' '.join(log_file9.split('_')).split('.')).split()
angle_scanned9 = name_split9[5]

log_file10 = sys.argv[10]
name_split10 = ' '.join(' '.join(log_file10.split('_')).split('.')).split()
angle_scanned10 = name_split10[5]

log_file11 = sys.argv[11]
name_split11 = ' '.join(' '.join(log_file11.split('_')).split('.')).split()
angle_scanned11 = name_split11[5]

log_file12 = sys.argv[12]
name_split12 = ' '.join(' '.join(log_file12.split('_')).split('.')).split()
angle_scanned12 = name_split12[5]

log_file13 = sys.argv[13]
name_split13 = ' '.join(' '.join(log_file13.split('_')).split('.')).split()
angle_scanned13 = name_split13[5]

log_file14 = sys.argv[14]
name_split14 = ' '.join(' '.join(log_file14.split('_')).split('.')).split()
angle_scanned14 = name_split14[5]

log_file15 = sys.argv[15]
name_split15 = ' '.join(' '.join(log_file15.split('_')).split('.')).split()
angle_scanned15 = name_split15[5]

log_file16 = sys.argv[16]
name_split16 = ' '.join(' '.join(log_file16.split('_')).split('.')).split()
angle_scanned16 = name_split16[5]

log_file17 = sys.argv[17]
name_split17 = ' '.join(' '.join(log_file17.split('_')).split('.')).split()
angle_scanned17 = name_split17[5]

log_file18 = sys.argv[18]
name_split18 = ' '.join(' '.join(log_file18.split('_')).split('.')).split()
angle_scanned18 = name_split18[5]


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
					cu_o_bond_dist1.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy1.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
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
				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0


with open(log_file6) as input:
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
					cu_o_bond_dist6.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy6.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0


with open(log_file7) as input:
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
					cu_o_bond_dist7.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy7.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0


with open(log_file8) as input:
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
					cu_o_bond_dist8.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy8.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0


with open(log_file9) as input:
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
					cu_o_bond_dist9.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy9.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0


with open(log_file10) as input:
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
					cu_o_bond_dist10.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy10.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0


with open(log_file11) as input:
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
					cu_o_bond_dist11.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy11.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0


with open(log_file12) as input:
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
					cu_o_bond_dist12.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy12.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0


with open(log_file13) as input:
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
					cu_o_bond_dist13.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy13.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0


with open(log_file14) as input:
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
					cu_o_bond_dist14.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy14.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0


with open(log_file15) as input:
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
					cu_o_bond_dist15.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy15.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0


with open(log_file16) as input:
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
					cu_o_bond_dist16.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy16.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0


with open(log_file17) as input:
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
					cu_o_bond_dist17.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy17.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0


with open(log_file18) as input:
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
					cu_o_bond_dist18.append(float(d[3]))
				if d[0]=="SCF" and d[1]=="Done:" and d[2]=="E(UB3LYP)" and d[3]=="=":
					#print d[4]
					QM_energy18.append(float(d[4]) * 627.503) 
				if d[0]=="Standard" and d[1]=="orientation:":
					count += 1
				if count==1:
					count2 += 1
				if count2 == n_atoms + 70:
					switch = 0
					count = 0
 					count2 = 0




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

first_element_in_qm_energy6 = QM_energy6[0]
for i in range(len(QM_energy6)):
	QM_energy6[i] -= first_element_in_qm_energy6

first_element_in_qm_energy7 = QM_energy7[0]
for i in range(len(QM_energy7)):
	QM_energy7[i] -= first_element_in_qm_energy7

first_element_in_qm_energy8 = QM_energy8[0]
for i in range(len(QM_energy8)):
	QM_energy8[i] -= first_element_in_qm_energy8

first_element_in_qm_energy9 = QM_energy9[0]
for i in range(len(QM_energy9)):
	QM_energy9[i] -= first_element_in_qm_energy9

first_element_in_qm_energy10 = QM_energy10[0]
for i in range(len(QM_energy10)):
	QM_energy10[i] -= first_element_in_qm_energy10

first_element_in_qm_energy11 = QM_energy11[0]
for i in range(len(QM_energy11)):
	QM_energy11[i] -= first_element_in_qm_energy11

first_element_in_qm_energy12 = QM_energy12[0]
for i in range(len(QM_energy12)):
	QM_energy12[i] -= first_element_in_qm_energy12

first_element_in_qm_energy13 = QM_energy13[0]
for i in range(len(QM_energy13)):
	QM_energy13[i] -= first_element_in_qm_energy13

first_element_in_qm_energy14 = QM_energy14[0]
for i in range(len(QM_energy14)):
	QM_energy14[i] -= first_element_in_qm_energy14

first_element_in_qm_energy15 = QM_energy15[0]
for i in range(len(QM_energy15)):
	QM_energy15[i] -= first_element_in_qm_energy15

first_element_in_qm_energy16 = QM_energy16[0]
for i in range(len(QM_energy16)):
	QM_energy16[i] -= first_element_in_qm_energy16

first_element_in_qm_energy17 = QM_energy17[0]
for i in range(len(QM_energy17)):
	QM_energy17[i] -= first_element_in_qm_energy17

first_element_in_qm_energy18 = QM_energy18[0]
for i in range(len(QM_energy18)):
	QM_energy18[i] -= first_element_in_qm_energy18


plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.plot(cu_o_bond_dist1[1:25], QM_energy1[0:24], "g^")
plt.plot(cu_o_bond_dist2[1:25], QM_energy2[0:24], "b^")
plt.plot(cu_o_bond_dist3[1:25], QM_energy1[0:24], "r^")
plt.plot(cu_o_bond_dist4[1:25], QM_energy1[0:24], "k^")
plt.plot(cu_o_bond_dist5[1:25], QM_energy1[0:24], "m^")
plt.plot(cu_o_bond_dist6[1:25], QM_energy1[0:24], "gs")
plt.plot(cu_o_bond_dist7[1:25], QM_energy1[0:24], "bs")
plt.plot(cu_o_bond_dist8[1:25], QM_energy1[0:24], "rs")
plt.plot(cu_o_bond_dist9[1:25], QM_energy1[0:24], "ks")
plt.plot(cu_o_bond_dist10[1:25], QM_energy1[0:24], "ms")
plt.plot(cu_o_bond_dist11[1:25], QM_energy1[0:24], "go")
plt.plot(cu_o_bond_dist12[1:25], QM_energy1[0:24], "bo")
plt.plot(cu_o_bond_dist13[1:25], QM_energy1[0:24], "ro")
plt.plot(cu_o_bond_dist14[1:25], QM_energy1[0:24], "ko")
plt.plot(cu_o_bond_dist15[1:25], QM_energy1[0:24], "mo")
plt.plot(cu_o_bond_dist16[1:25], QM_energy1[0:24], "g+")
plt.plot(cu_o_bond_dist17[1:25], QM_energy1[0:24], "b+")
plt.plot(cu_o_bond_dist18[1:25], QM_energy1[0:24], "r+")

plt.legend([angle_scanned1, angle_scanned2, angle_scanned3, angle_scanned4, angle_scanned5, angle_scanned6, angle_scanned7, angle_scanned8, angle_scanned9, angle_scanned10, angle_scanned11, angle_scanned12, angle_scanned13, angle_scanned14, angle_scanned15, angle_scanned16, angle_scanned17, angle_scanned18], fontsize='10', bbox_to_anchor=(.75, .9, 0.25, 0.0), loc=3, ncol=1, mode='expand', borderaxespad=0., numpoints = 1)


plt.ylabel('Energy (kcal/mol)')
plt.xlabel('Copper-Oxygen Distance ($\AA$)')
plt.xlim((0,4))
plt.axhline(y = 0, color = "k")
plt.title(r'Comparing B3LYP and PBE1PBE functionals of the 5-Water System', size='14')
#plt.suptitle(r'With Parameters from 1 Water System', size='14')
plt.show()







