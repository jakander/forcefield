import sys
import math
import numpy as np
import matplotlib.pyplot as plt


bond_dist = [1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0]
sigma = 1.65
epsilon = 17.4
LJ_energy = []
for r in bond_dist:
	LJ= 4*epsilon*((sigma/r)**12 - (sigma/r)**6)
	LJ_energy.append(float(LJ))

#compute the sigma and epsilon values for the Lennard Jones through a least squares approach
repulsive_term = [12]
colors_repulsive_term = ['bs']
for m in range(len(repulsive_term)):
	r1 = []
	r2 = []
	r3 = []
	for r in bond_dist:
		l = r ** (-repulsive_term[m])
		k = -r ** (-6)
		r1.append(l)
		r2.append(k)
		r3.append(1)
	r_all = [r1, r2, r3]
	r_all = np.asmatrix(r_all)
	LJ_energy = np.asarray(LJ_energy)
	A, B, C = np.linalg.lstsq(r_all.T, LJ_energy)[0]
	residuals = np.linalg.lstsq(r_all.T,LJ_energy)[1]
	print A, B, C
	print"for m =", (m+12),"sum of residuals =", residuals
	math = .3333333333333333333333333333333333333333333/2 
	sigma1 = (A/B)**math
	epsilon1 = (B**2) / (4*A)
	print "sigma =", sigma1, ",  ", "epsilon =", epsilon1
	Fit_LJ = []
	for r in bond_dist:
		l = 4*epsilon1*((sigma1/r)**12 - (sigma1/r)**6) 
		Fit_LJ.append(l)
	#element1inFitLJ = Fit_LJ[0]
	#for i in range(len(Fit_LJ)):
	#	Fit_LJ[i] -= element1inFitLJ 
	plt.plot(bond_dist, Fit_LJ, colors_repulsive_term[m])


plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.plot(bond_dist[1:26], LJ_energy[0:25], "r^")
plt.ylabel('Energy (kcal/mol)')
plt.xlabel('Copper-Oxygen Distance ($\AA$)')
plt.legend(['LJ_fitted', 'LJ'], fontsize='10', bbox_to_anchor=(.84, 0.905, 0.16, 0.0), loc=3, ncol=1, mode='expand', borderaxespad=0., numpoints = 1)
plt.xlim((0,4))
plt.axhline(y = 0, color = "k")
plt.title(r'Testing LJ fitting', size='14')
#plt.suptitle(r'With Parameters from 1 Water System', size='14')
plt.show()
