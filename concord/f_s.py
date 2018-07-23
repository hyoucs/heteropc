import numpy as np
from scipy.io import loadmat
from scipy.io import savemat
from numpy.linalg import matrix_power
from numpy.linalg import inv

from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--mode", "-m", type=int, help="Mode", default=2)
parser.add_argument("--pwr", "-p", type=int, help="Power", default=10)
parser.add_argument("--pmode", "-pm", type=bool, help="Power Mode", default=False)
args = parser.parse_args()

def func1(x, a, b, c):
	return a * np.exp(-b * x) + c
def func2(x, a, b, c):
	return a*x/(b+x)+c
def cookD(x, e):
	x = np.asarray(x).reshape(-1,1)
	num = x.shape[0]
	t = x - np.mean(x)
	h = t**2/sum(sum(t**2))+1/num

	mse = sum(e**2)/(num-1)
	# d = np.sign(e)*h*e**2/(3*mse*(1-h)**2)
	d = np.zeros(num)
	for i in range(num):
		d[i] = np.sign(e[i])*h[i]*e[i]**2/(3*mse*(1-h[i])**2)
	return d

Dir = '/Users/LSK/Dropbox/glasso/data/HCP-V1/'

sMat = loadmat(Dir+'Diffusion-q.mat')
s = sMat['X']
fMat = loadmat(Dir+'tfMRI-EMOTION.mat')
f = fMat['X']

d,_,n = s.shape

mode = args.mode
pwr = args.pwr
pwrMode = args.pmode
pwrS = np.zeros((d,d,n))
for i in range(n):
	pwrS[:,:,i] = matrix_power(s[:,:,i],pwr)
# check (i,j)th entry across 51 subjects
if mode == 1:
	x = []
	y = []
	z = []
	ctr = 0
	for i in range(d-1):
		for j in range(i+1,d):
			for k in range(n):
				x.append(s[i,j,k])
				y.append(f[i,j,k])
				z.append(pwrS[i,j,k])
			x = np.asarray(x)
			y = np.asarray(y)
			z = np.asarray(z)

			popt, _ = curve_fit(func2, x, y)
			plt.plot(x,y,'bo')
			plt.plot(x, func2(x, *popt), 'ro')
			if pwrMode:
				plt.plot(z,y,'go')
			plt.ylim(-1.0,1.0)
			plt.show()
			x = []
			y = []
			z = []
			ctr = ctr+1
			if ctr==10:
				exit()

# check all entries of the same subject
if mode == 2:
	x = []
	y = []
	z = []
	# ctr = 0
	L = []

	for k in range(n):
		for i in range(d-1):
			for j in range(i+1,d):
				# if(s[i,j,k]>100):
				x.append(s[i,j,k])
				y.append(f[i,j,k])
				z.append(pwrS[i,j,k])
		x = np.asarray(x)
		y = np.asarray(y)
		z = np.asarray(z)

		popt, _ = curve_fit(func2, z, y)
		e = y - func2(z, *popt)

		# r = (abs(cookD(x,e))*10000).tolist()
		# plt.scatter(x,y,s=r)
		# plt.plot(x, func2(x, *popt), 'ro')
		# if pwrMode:
		# 	plt.plot(z,y,'go')
		# plt.show()
		# L.append(cookD(x,e))
		L.append(e)
		x = []
		y = []
		z = []
		# ctr = ctr+1
		# if ctr==10:
		# 	exit()
	r = np.vstack(L)
	a = {}
	a['X'] = r
	savemat('sf_e.mat',a)