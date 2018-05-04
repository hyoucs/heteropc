"""Drawing the overlapped network across different functions"""
import numpy as np
from scipy.io import loadmat

import matplotlib.pyplot as plt
import networkx as nx

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def plot(omega):
	r_name = [
	'lateralorbitofrontal-l',
	'parsorbitalis-l',
	'frontalpole-l',
	'medialorbitofrontal-l',
	'parstriangularis-l',
	'parsopercularis-l',
	'rostralmiddlefrontal-l',
	'superiorfrontal-l',
	'caudalmiddlefrontal-l',
	'precentral-l',
	'paracentral-l',
	'rostralanteriorcingulate-l',
	'caudalanteriorcingulate-l',
	'posteriorcingulate-l',
	'isthmuscingulate-l',
	'postcentral-l',
	'supramarginal-l',
	'superiorparietal-l',
	'inferiorparietal-l',
	'precuneus-l',
	'cuneus-l',
	'pericalcarine-l',
	'lateraloccipital-l',
	'lingual-l',
	'fusiform-l',
	'parahippocampal-l',
	'entorhinal-l',
	'temporalpole-l',
	'inferiortemporal-l',
	'middletemporal-l',
	'bankssts-l',
	'superiortemporal-l',
	'transversetemporal-l',
	'insula-l',
	'thalamusproper-l',
	'caudate-l',
	'putamen-l',
	'pallidum-l',
	'accumbensarea-l',
	'hyppocampus-l',
	'amygdala-l',
	'lateralorbitofrontal-r',
	'parsorbitalis-r',
	'frontalpole-r',
	'medialorbitofrontal-r',
	'parstriangularis-r',
	'parsopercularis-r',
	'rostralmiddlefrontal-r',
	'superiorfrontal-r',
	'caudalmiddlefrontal-r',
	'precentral-r',
	'paracentral-r',
	'rostralanteriorcingulate-r',
	'caudalanteriorcingulate-r',
	'posteriorcingulate-r',
	'isthmuscingulate-r',
	'postcentral-r',
	'supramarginal-r',
	'superiorparietal-r',
	'inferiorparietal-r',
	'precuneus-r',
	'cuneus-r',
	'pericalcarine-r',
	'lateraloccipital-r',
	'lingual-r',
	'fusiform-r',
	'parahippocampal-r',
	'entorhinal-r',
	'temporalpole-r',
	'inferiortemporal-r',
	'middletemporal-r',
	'bankssts-r',
	'superiortemporal-r',
	'transversetemporal-r',
	'insula-r',
	'thalamusproper-r',
	'caudate-r',
	'putamen-r',
	'pallidum-r',
	'accumbensarea-r',
	'hyppocampus-r',
	'amygdala-r',
	'brainstem'
	]

	[d,_] = omega.shape #dim: 3403*3403
	r = len(r_name) #num_regions

	## Define dict
	idx = {}
	ctr = 0
	for i in range(r-1):
		for j in range(i+1,r):
			idx[i,j] = ctr
			ctr = ctr + 1

	## Define inverse dict
	inv_idx = {v: k for k, v in idx.iteritems()}

	G = nx.Graph()
	G.add_nodes_from(np.arange(r))

	ctr = 0
	for i in range(d-1):
		for j in range(i+1,d):
			if omega[i,j] != 0:
				G.add_edge(inv_idx[i][0],inv_idx[i][1])
				G.add_edge(inv_idx[j][0],inv_idx[j][1])
				ctr = ctr+1
	print(ctr)
	## Delete nodes with degree 0
	to_del = nx.isolates(G)
	G.remove_nodes_from(to_del)

	## Draw
	options = {'node_color': '#FA8072', 'edge_color': '#2C3E50', 
				'node_size': 400,'width': 0.8,}
	labels = {}
	for i in range(r):
		labels[i]=r_name[i]
	G = nx.relabel_nodes(G,labels)
	nx.draw(G, with_labels=True, **options) #font_weight='bold'
	plt.show()

##### Plot the common network from 2 function results ####
mat1 = loadmat('/Users/LSK/Dropbox/glasso/concord_results/GAMBLING/c/0.28__s.mat')
invcov1 = mat1['S']
invcov1 = invcov1.toarray()
mat2 = loadmat('/Users/LSK/Dropbox/glasso/concord_results/EMOTION/no_c/0.295_s.mat')
invcov2 = mat2['S']
invcov2 = invcov2.toarray()

[p,_] = invcov1.shape
omega1 = np.zeros((p,p))
for i in range(p):
	for j in range(p):
		if(invcov1[i,j] != 0 and invcov2[i,j] != 0):
			omega1[i,j] = 1
plot(omega1)

##### Plot network from a single result #####
mat = loadmat('/Users/LSK/Dropbox/glasso/concord_results/cross_func/g_e(114)/0.311_114.mat')
omega2 = mat['S']
omega2 = omega2.toarray()
plot(omega2)

##### Plot the common network found by both methods #####
omega = np.zeros((p,p))
for i in range(p):
	for j in range(p):
		if(omega1[i,j] != 0 and omega2[i,j] != 0):
			omega[i,j] = 1
plot(omega)