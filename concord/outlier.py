import numpy as np
import math
from scipy.io import loadmat
import os

def check_symmetric(a, tol=1e-8):
    return np.allclose(a, a.T, atol=tol)

avg_list = []
dir = '/Users/LSK/Dropbox/glasso/concord_results/cross_func/l_r/'
for filename in os.listdir(dir):
	if filename.endswith(".mat"): 
		print(filename)

		f_n = os.path.join(dir, filename)
		mat = loadmat(f_n)
		omega = mat['S']
		omega = omega.toarray()

		# Partial correlation
		rho = np.zeros(omega.shape)
		d = omega.shape[0]

		for i in range(d):
			for j in range(i,d):
				rho[i,j] = -omega[i,j]/math.sqrt(omega[i,i]*omega[j,j])
		i_lower = np.tril_indices(d, -1)
		rho[i_lower] = rho.T[i_lower]
		# print(check_symmetric(rho))

		np.fill_diagonal(rho, 0.0)
		nz = np.count_nonzero(rho)
		print(nz)

		avg = sum(sum(np.absolute(rho)))/nz
		avg_list.append(avg)
		print(avg)

print(sum(avg_list)/float(len(avg_list)))