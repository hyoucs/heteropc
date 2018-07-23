import numpy as np
import numpy.random as rng
from scipy.io import loadmat
from numpy.linalg import inv, norm
import math
import cvxpy as cvx

######### Data Preparation #########
# Dir = '/Users/LSK/Dropbox/glasso/data/HCP-V1/'
# sMat = loadmat(Dir+'Diffusion-q.mat')
# s = sMat['X']
# fMat = loadmat(Dir+'tfMRI-EMOTION.mat')
# f = fMat['X']

# d,_,n = s.shape

# vec_s = []
# vec_f = []
# p = 0
# for i in range(d-1):
# 	for j in range(i+1,d):
# 		vec_s.append(s[i,j])
# 		vec_f.append(f[i,j])
# 		p = p+1
# vec_s = np.transpose(np.asarray(vec_s))
# vec_f = np.transpose(np.asarray(vec_f))
####################################
Dir = '/Users/LSK/Dropbox/glasso/heteropc/concord/synthetic.mat'
mat_data = loadmat(Dir)
vec_s = mat_data['S']
vec_f = mat_data['F']
n,p = vec_s.shape
####################################
print(vec_s.shape)
print(vec_f.shape)
# Missing entry indices
r_miss = 0.2
# miss = rng.choice(p,int(r_miss*p),replace=False)
# known = [x for x in np.arange(p) if x not in miss
q = int(r_miss*p) # miss number
K = 10
test_num = int(n/K)

######### Original model error #########
err1 = []
for i in range(K):
	# Load model i's Omega gotten from R script
	dir_i = 'syn/f/'+str(i)+'.mat'
	mat = loadmat(dir_i)
	Omega = mat['S']
	# Omega = Omega.toarray()

	testD = vec_f[i*test_num:(i+1)*test_num,:]
	d1 = vec_f[:i*test_num]
	d2 = vec_f[(i+1)*test_num:]
	trainD = np.vstack((d1,d2))
	# Model i's mu vector
	mu = np.mean(trainD, axis=0)
	Sigma = inv(Omega)
	# completion given partial correlation
	mu1 = mu[:q]
	mu2 = mu[q:]
	S12 = Sigma[:q,q:]
	Omega22 = Omega[q:,q:]

	err = 0
	for j in range(test_num):
		a = testD[j,q:]
		_mu = mu1+np.matmul(S12,np.matmul(Omega22,a-mu2))
		err = err + math.sqrt(sum((_mu-testD[j,:q])**2))/math.sqrt(sum(testD[j,:q]**2))
	# print err/test_num
	err1.append(err/test_num)
err1 = np.mean(err1)
print(err1)

######### Regression bsased completion error #########
# err2 = []
# for k in range(K):
# 	# Load model k's Omega gotten from R script
# 	dir_i = '10fold/e/'+str(k)+'.mat'
# 	mat = loadmat(dir_i)
# 	Omega = mat['S']
# 	Omega = Omega.toarray()

# 	diag = np.repeat(Omega.diagonal().reshape((p,1)),p,axis=1)
# 	beta = -np.divide(Omega,diag)

# 	testD = vec_f[k*test_num:(k+1)*test_num,:]
	
# 	# now do regression one test sample by one
# 	err = 0
# 	for w in range(test_num):
# 		y_a = cvx.Variable(q)
# 		y_c = testD[w,q:]

# 		loss = 0
# 		for i in range(p):
# 			if i in range(q):
# 				loss_ = y_a[i]
# 				for j in range(p-q):
# 					loss_ = loss_-beta[i,q+j]*y_c[j]
# 				for j in range(q):
# 					if (j != i):
# 						loss_ = loss_-beta[i,j]*y_a[j]
# 			else:
# 				loss_ = y_c[i-q]
# 				for j in range(q):
# 					loss_ = loss_-beta[i,j]*y_a[j]
# 				for j in range(p-q):
# 					if (j+q != i):
# 						loss_ = loss_-beta[i,j+q]*y_c[j]
# 			loss_n = loss_**2
# 			loss = loss+loss_n
# 		problem = cvx.Problem(cvx.Minimize(loss))
# 		problem.solve()

# 		err = err + math.sqrt(sum((y_a-testD[w,:q])**2))/math.sqrt(sum(testD[w,:q]**2))
# 	print err/test_num
# 	err2.append(err/test_num)
# err2 = np.mean(err2)
# print(err2)

######### S_F model error #########
err3 = []
vec_s = vec_s/np.amax(vec_s)
for i in range(K):
	dir_i = 'syn/sf/'+str(i)+'.mat'
	mat = loadmat(dir_i)
	Omega = mat['S']
	# Omega = Omega.toarray()

	test_f = vec_f[i*test_num:(i+1)*test_num,:]
	test_s = vec_s[i*test_num:(i+1)*test_num,:]
	testD = test_f-test_s

	train_f = np.vstack((vec_f[:i*test_num],vec_f[(i+1)*test_num:]))
	train_s = np.vstack((vec_s[:i*test_num],vec_s[(i+1)*test_num:]))
	# trainD = np.sign(train_f)*(abs(train_f)-abs(train_s))
	trainD = train_f-train_s

	mu = np.mean(trainD, axis=0)
	Sigma = inv(Omega)

	mu1 = mu[:q]
	mu2 = mu[q:]
	S12 = Sigma[:q,q:]
	Omega22 = Omega[q:,q:]

	err = 0
	for j in range(test_num):
		a = testD[j,q:]
		_mu = mu1+np.matmul(S12,np.matmul(Omega22,a-mu2))
		err = err + math.sqrt(sum((_mu+test_s[j,:q]-test_f[j,:q])**2))/math.sqrt(sum(test_f[j,:q]**2))
	# print err/test_num
	err3.append(err/test_num)
err3 = np.mean(err3)
print(err3)