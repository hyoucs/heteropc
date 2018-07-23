import numpy as np
from scipy.io import savemat, loadmat

root = 'result/'

def create_pMat(arg,p=210):
	assert arg=='sf_ss' or 'ff' or 'ss' or 'all', 'Invalid argument'
	q = int(p/2)

	pMat = np.zeros((p,p))
	np.fill_diagonal(pMat,1)

	if arg == 'sf_ss':
		pMat[:,q:] = 1
		pMat[q:,:] = 1
	elif arg == 'ff':
		pMat[:q,:q] = 1
	elif arg == 's_s':
		pMat[q:,q:]=1
	elif arg == 'all':
		pMat[:,:] = 1

	a = {}
	a['M'] = pMat
	fname = 'pMask/pMask_'+arg+'.mat'
	savemat(fname,a)

def check_symmetric(a, tol=1e-8):
	return np.allclose(a, a.T, atol=tol)

def test_edge(omega,diff=False):
	assert check_symmetric(omega)==True, 'Input must be symmetric'
	p,_ = omega.shape
	if diff:
		print(int(np.count_nonzero(omega)/2))
	else:
		print(int((np.count_nonzero(omega)-p)/2))
	F_F = 0
	S_F = 0
	S_S = 0
	q = 105
	for i in range(p-1):
		for j in range(i+1,p):
			if omega[i,j] != 0:
				if j <= q:
					F_F += 1
				elif i <= q:
					S_F += 1
				else:
					S_S += 1
	print('F_F has {} edges, S_F has {} edges, S_S has {} edges'.format(F_F,S_F,S_S))

def test_diff(dir1,dir2):
	matdata = loadmat(dir1)
	omega1 = matdata['M']
	omega1[omega1!=0] = 1
	matdata = loadmat(dir2)
	omega2 = matdata['M']
	omega2[omega2!=0] = 1
	test_edge(omega1-omega2,diff=True)

def test_one(dir_):
	matdata = loadmat(dir_)
	omega = matdata['M']
	test_edge(omega)

if __name__ == '__main__':
	arg = 'ss'
	create_pMat(arg)
	# test_one(root+'0.259_702.mat')
	# test_diff(root+'0.259_702.mat',root+'0.259_580.mat')
	