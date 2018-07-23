import numpy as np
import numpy.random as rng
from scipy.io import savemat
from numpy.linalg import inv,eigh
from scipy.stats import truncnorm

from gen_from_cor import gen_from_cor

def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
	return truncnorm(
		(low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)
def nz(a):
	return np.sum(a!=0)
def check_symmetric(a, tol=1e-8):
	return np.allclose(a, a.T, atol=tol)
def isPSD(A, tol=1e-8):
	E,V = eigh(A)
	return np.all(E > -tol)

# some constants
p = 105
n = 500
nz_ratio = 0.1
nz_val = 0.8
rmean = 0.6
rsd = 0.3

# covariance matrix from precision matrix
d = p*(p-1)/2
num_nz = int(nz_ratio*d)
idx = rng.choice(d, num_nz, replace=False)
u = np.zeros(d)
u[idx] = nz_val

Omega = np.zeros((p,p))
Omega[np.triu_indices(p,1)] = u
Omega = Omega+Omega.T
np.fill_diagonal(Omega,1)
Sigma = inv(Omega)

# generate 'function': shape n*p
mean = np.zeros(p)
f = rng.multivariate_normal(mean,Sigma,n)
print(f.shape)

# correlations
x = get_truncated_normal(mean=rmean, sd=rsd, low=0, upp=1)
r = x.rvs(p)
# randomly flip 50% r into negative values
# idx = rng.choice(p,int(0.5*p),replace=False)
# r[idx] = -r[idx]

# generate 'structures'
slist = []
for i in range(p):
	slist.append(gen_from_cor(f[:,i],r[i]).reshape(-1,1))
s = np.hstack(slist)
print(s.shape)

a = {}
a['S'] = s
a['F'] = f
savemat('synthetic.mat',a)