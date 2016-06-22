import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

stdev = np.std
sqrt = np.sqrt
nullfmt = NullFormatter()


dat1 = sys.argv[1]
system = sys.argv[2]


def hist1d(data, x_axis, num_b, system, analysis, norm):
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
	
	events, edges, patches = plt.hist(data, bins=num_b, histtype = 'bar', normed=norm)
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel(r'%s' %(x_axis), size=12)
#	plt.title(r'something here...' %(system), size='14')

	if norm == True:
		plt.ylabel('Probability Density')
		plt.savefig('%s.%s.prob1d.png' %(system,analysis))
		nf = open('%s.%s.prob1d.dat' %(system,analysis),'w')
	else:
		plt.ylabel('Frequency', size=12)
		plt.savefig('%s.%s.hist1d.png' %(system,analysis))
		nf = open('%s.%s.hist1d.dat' %(system,analysis), 'w')

	for i in range(len(events)):
		nf.write('%10.1f      %10.4f\n' %(events[i], edges[i]))
	
	plt.close()
	nf.close()
	events = []
	edges = []
	patches = []













