#Checks if a particular Xij has a normal distribution over all 51 subjects

import os
import numpy as np
import scipy
import scipy.stats
import scipy.io as spio
import codecs

#All distributions in scipy
dist_names = [ 'alpha', 'anglit', 'arcsine', 'beta', 'betaprime', 'bradford', 'burr', 'cauchy', 'chi', 'chi2', 'cosine', 'dgamma', 'dweibull', 'erlang', 'expon', 'exponweib', 'exponpow', 'f', 'fatiguelife', 'fisk', 'foldcauchy', 'foldnorm', 'frechet_r', 'frechet_l', 'genlogistic', 'genpareto', 'genexpon', 'genextreme', 'gausshyper', 'gamma', 'gengamma', 'genhalflogistic', 'gilbrat', 'gompertz', 'gumbel_r', 'gumbel_l', 'halfcauchy', 'halflogistic', 'halfnorm', 'hypsecant', 'invgamma', 'invgauss', 'invweibull', 'johnsonsb', 'johnsonsu', 'ksone', 'kstwobign', 'laplace', 'logistic', 'loggamma', 'loglaplace', 'lognorm', 'lomax', 'maxwell', 'mielke', 'nakagami', 'ncx2', 'ncf', 'nct', 'norm', 'pareto', 'pearson3', 'powerlaw', 'powerlognorm', 'powernorm', 'rdist', 'reciprocal', 'rayleigh', 'rice', 'recipinvgauss', 'semicircular', 't', 'triang', 'truncexpon', 'truncnorm', 'tukeylambda', 'uniform', 'vonmises', 'wald', 'weibull_min', 'weibull_max', 'wrapcauchy']

filenames = [each.rstrip('\n') for each in os.listdir() if each.endswith(".mat")]
counter = np.zeros((83,83))

for f in filenames:
    try:
        mat=spio.loadmat(f, squeeze_me=True)
        print("Loading file ", f)
    except:
        print("Bad file")
        continue
    for i in range(83):
        for j in range(i+1.83):
            #print("Xij value for i=",i,"j=",j)
            x1, p = (scipy.stats.normaltest(mat['X'][i][j]))
            if (p>0.04):
                counter[i][j] += 1
                #print("Xij value for i=",i,"j=",j)
                #print(x1,p)

print(counter)
