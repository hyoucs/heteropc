"""Drawing the overlapped network across different functions"""
import numpy as np
from scipy.io import loadmat

import matplotlib.pyplot as plt
import networkx as nx

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

mat1 = loadmat('/Users/LSK/Dropbox/glasso/concord_results/EMOTION/emotion.mat')
invcov1 = mat1['invcov']
mat2 = loadmat('/Users/LSK/Dropbox/glasso/concord_results/RELATIONAL/relational.mat')
invcov2 = mat2['invcov']

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
nx.draw(G, with_labels=True, **options) #font_weight='bold'
plt.show()