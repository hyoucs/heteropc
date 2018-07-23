import numpy as np
from math import acos,tan
from numpy.linalg import qr
import numpy.random as rng

def gen_from_cor(x1,rho):
	theta = acos(rho)
	x1 = x1.reshape(-1,1)
	n = x1.shape[0]
	x2 = (0.5*rng.randn(n)+0.5).reshape(-1,1)
	X = np.hstack((x1,x2))
	Xctr = X-np.mean(X,axis=0)

	Id = np.diag(np.ones(n))
	Q,_ = qr(Xctr[:,0].reshape(-1,1))

	P = Q.dot(Q.T)
	x2o = (Id-P).dot(Xctr[:, 1]).reshape(-1,1)
	Xc2 = np.hstack((Xctr[:, 0].reshape(-1,1), x2o))
	Y = Xc2.dot(np.diag(1.0/np.sqrt(np.sum(Xc2**2,axis=0))))
	x = Y[:,1]+(1/tan(theta))*Y[:,0]
	# print(np.corrcoef(x1.flatten(),x))
	return x

if __name__ == '__main__':
	n = 20
	rho = 0.6
	x1 = rng.randn(n)+1
	x = gen_from_cor(x1,rho)
	print(np.corrcoef(x1,x))