
#PREAMBLE:

import numpy as np
from scipy import stats
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

stdev = np.std
sqrt = np.sqrt
nullfmt = NullFormatter()


dat1 = sys.argv[1]
system = sys.argv[2]
sel = []
sel.append(['bond dist','firebrick'])

t = "Mean= 2.044"\
	" Mode= 2.016"


#-----------------------------------------------------------------------
#Subroutine:


def hist1d(data, x_axis, num_b, system, norm):
	""" Creates a 1D histogram:
	Usage: hist1d(data, x_axis, num_b, system, analysis, norm)
	
	Arguments:
	data: self-explanatory
	x_axis: string to be used for the axis label
	num_b: number of bins to be used when binning the data
	system: descriptor for the system analyzed
	analysis: descriptor for the analysis performed and plotted
	norm = [False][True]; if False, plotting a frequency of data; if True, plotting a probability density
	"""
	mode = stats.mode(data)
	mean = np.mean(data)
	stdev = np.std(data)
	print 'stdev=', stdev
	
	plt.hist(data, bins=num_b, histtype = 'bar', normed=norm)
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel(r'%s' %(x_axis), size=12)
	plt.title('From CSD for All Structures with Cu Coordinated to Two N', size='14')
	plt.axvline(x = mean, color = "r")
	plt.xticks(np.arange(1.3, 3.5, 0.2))
	plt.xlim((1.7,2.7))
	#plt.text(2.8, 8000, t, ha='left', va='top', wrap=True, fontsize=12, color='k', bbox=dict(facecolor='#D3D3D3', ec='none'))
	#plt.text(3.0,7000, "Mode=2.016", va='top', ha='right', fontsize=12, color='k', bbox=dict(facecolor='#D3D3D3',ec='none'))
	plt.text(2.5,8500, "Mean= 2.044" , va='top', ha='left', fontsize=12, color='k', backgroundcolor='#D3D3D3')
	plt.text(2.5,8000, "Mode= 2.016", va='top', ha='left', fontsize=12, color='k', backgroundcolor='#D3D3D3')

	if norm == True:
		plt.ylabel('Probability Density')
		#plt.savefig('%s.%s.prob1d.png' %(system))
		#nf = open('%s.%s.prob1d.dat' %(system),'w')
	else:
		plt.ylabel('Frequency', size=12)
		#plt.savefig('%s.%s.hist1d.png' %(system))
		#nf = open('%s.%s.hist1d.dat' %(system), 'w')
	
	
	plt.show()
	plt.close()




#-------------------------------------------------------------------------
#MAINPROGRAM:


datalist1 = np.loadtxt(dat1)
data_range = np.percentile(datalist1[:,:], [0, 100])
print data_range


for i in range(len(sel)):
	selection = sel[i][0]

	hist1d(datalist1[:,i], 'Coordination Distance ($\AA$)', 250, '%02d.%s.%s' %(i, system, selection), False)













