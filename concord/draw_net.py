"""Drawing the overlapped network across different functions"""
import numpy as np
from scipy.io import loadmat

import matplotlib.pyplot as plt
import networkx as nx

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

mat1 = loadmat('/Users/LSK/Dropbox/glasso/concord_results/RELATIONAL/relational.mat')
invcov1 = mat1['invcov']
mat2 = loadmat('/Users/LSK/Dropbox/glasso/concord_results/GAMBLING/gambling.mat')
invcov2 = mat2['invcov']

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

[d,_] = invcov1.shape #dim: 3403*3403
r = 83 #num_regions

## Define dict
idx = {}
ctr = 0
for i in range(r-1):
	for j in range(i+1,r):
		idx[i,j] = ctr
		ctr = ctr + 1

## Define inverse dict
inv_idx = {v: k for k, v in idx.iteritems()}

# invcov1[np.where(invcov1!=0)]=1
# invcov2[np.where(invcov2!=0)]=1
# print sum(sum(abs(invcov2-invcov1)))
G = nx.Graph()
G.add_nodes_from(np.arange(r))

for i in range(d-1):
	for j in range(i+1,d):
		if(invcov1[i,j] != 0 and invcov2[i,j] != 0):
			G.add_edge(inv_idx[i][0],inv_idx[i][1])
			G.add_edge(inv_idx[j][0],inv_idx[j][1])
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