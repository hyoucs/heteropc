import numpy as np 
import numpy.random as rng
from scipy.io import savemat

# Common precision matrix
# Build a dict for other entries
p = 20
val = 0.3
base_size = 3
diff_size = 20
diff_val = 0.5

base = np.identity(p)
idx = {}
ctr = 0
for i in range(p):
	for j in range(i+1,p):
		if i in range(base_size) and j in range(p-base_size,p):
			base[i,j] = val
			base[j,i] = val
		else:
			idx[ctr] = [i,j]
			ctr = ctr+1

z = len(idx)
print(z)

# create 2 omega using the same base
o_idx = rng.choice(z,diff_size*2,replace=False)

omega1 = base
o1_idx = o_idx[:diff_size]
for i in range(len(o1_idx)):
	cur_idx = idx[o1_idx[i]]
	omega1[cur_idx[0],cur_idx[1]] = diff_val
	omega1[cur_idx[1],cur_idx[0]] = diff_val

omega2 = base
o2_idx = o_idx[-diff_size:]
for i in range(len(o2_idx)):
	cur_idx = idx[o2_idx[i]]
	omega2[cur_idx[0],cur_idx[1]] = diff_val
	omega2[cur_idx[1],cur_idx[0]] = diff_val

# sample data from omega1 and omega2
mean = np.zeros(p)

sigma1 = np.linalg.inv(omega1)
data1 = rng.multivariate_normal(mean,sigma1,500)
a = {}
a['M'] = data1
savemat('outlier/data_sw1.mat', a) 
sigma2 = np.linalg.inv(omega2)
data2 = rng.multivariate_normal(mean,sigma2,500)
a = {}
a['M'] = data2
savemat('outlier/data_sw2.mat', a)