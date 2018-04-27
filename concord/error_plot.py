''' 
Plot error (likelihood and L1 penalty) vs number of edge pair/ lambda curve
'''
import numpy as np
from scipy.io import loadmat

import matplotlib.pyplot as plt

# Const to change
lam = 0.8

# Data
mat1 = loadmat('/Users/LSK/Dropbox/glasso/concord_results/ASD/0.8_14_s.mat')
Omega = mat1['S']
Omega = Omega.toarray()
ctr = 0
pos = 0
for i in range(np.shape(Omega)[0]):
	for j in range(np.shape(Omega)[1]):
		if Omega[i,j]!=0 and i!=j:
			ctr = ctr+1
			if Omega[i,j]>0:
				pos = pos+1
			# print(str(i)+','+str(j)+': '+str(Omega[i,j]))
print(Omega[:5,:5])
print(pos)
print(ctr)

pMask = loadmat('/Users/LSK/Dropbox/glasso/pMask.mat')
pMat = pMask['M'] #3403*3403
pMat[pMat==0] = 100
np.fill_diagonal(pMat, 0)

mat2 = loadmat('/Users/LSK/Dropbox/glasso/concord_results/EMOTION/data.mat')
D = mat2['M']
# Standardize
S = D-np.tile(np.mean(D, axis=0),(np.shape(D)[0],1))
S = np.matmul(np.transpose(S),S) / (np.shape(S)[0] - 1)

# negative Likelihoods
term1 = -sum(np.log(Omega.diagonal()))
term2 = 0.5*np.trace(np.matmul(Omega,np.matmul(S,Omega)))
print(term1)
print(term2)
print(term1+term2)

# L1 penalty
p_term = np.multiply(np.abs(Omega),lam*pMat)
L1_with_lam = sum(sum(p_term))
L1_penalty = sum(sum(np.abs(Omega)))
print(L1_with_lam)
print(L1_penalty)

print(lam)