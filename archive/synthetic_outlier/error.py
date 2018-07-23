''' 
Plot error (likelihood and L1 penalty) vs number of edge pair/ lambda curve
'''
import numpy as np
from scipy.io import loadmat

# Const to change
lam = 0.2
# Data
mat1 = loadmat('c/0.2_138.mat')
Omega = mat1['M']
ctr = 0
pos = 0
for i in range(np.shape(Omega)[0]):
	for j in range(np.shape(Omega)[1]):
		if Omega[i,j]!=0 and i!=j:
			ctr = ctr+1
			if Omega[i,j]>0:
				pos = pos+1
print(Omega[:5,:5])
print(pos)
print(ctr)

mat2 = loadmat('lw/data_lw.mat')
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
p_term = np.abs(Omega)*lam
L1_with_lam = sum(sum(p_term))
L1_penalty = sum(sum(np.abs(Omega)))
print(L1_with_lam)
print(L1_penalty)

print(lam)