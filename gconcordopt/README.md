# README #

### Installation ###

Rcpp and RcppEigen may need to be installed separately

```
#!bash

syoh@skipper:~/Downloads$ hg clone ssh://hg@bitbucket.org/sangoh/gconcordopt
remote: Warning: Permanently added the RSA host key for IP address '104.192.143.3' to the list of known hosts.
destination directory: gconcordopt
requesting all changes
adding changesets
adding manifests
adding file changes
added 30 changesets with 100 changes to 29 files (+1 heads)
updating to branch default
12 files updated, 0 files merged, 0 files removed, 0 files unresolved

syoh@skipper:~/Downloads$ R CMD build gconcordopt 
* checking for file ‘gconcordopt/DESCRIPTION’ ... OK
* preparing ‘gconcordopt’:
* checking DESCRIPTION meta-information ... OK
* cleaning src
* checking for LF line-endings in source and make files
* checking for empty or unneeded directories
* building ‘gconcordopt_1.0.tar.gz’

syoh@skipper:~/tmp$ R CMD INSTALL gconcordopt_1.0.tar.gz 
* installing to library ‘/home/syoh/R/x86_64-pc-linux-gnu-library/3.2’
* installing *source* package ‘gconcordopt’ ...
** libs
g++ -I/usr/share/R/include -DNDEBUG -DEIGEN_NO_DEBUG  -I"/usr/local/lib/R/site-library/Rcpp/include" -I"/usr/local/lib/R/site-library/RcppEigen/include"   -fpic  -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c RcppExports.cpp -o RcppExports.o
g++ -I/usr/share/R/include -DNDEBUG -DEIGEN_NO_DEBUG  -I"/usr/local/lib/R/site-library/Rcpp/include" -I"/usr/local/lib/R/site-library/RcppEigen/include"   -fpic  -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c gconcordopt.cpp -o gconcordopt.o
g++ -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions -Wl,-z,relro -o gconcordopt.so RcppExports.o gconcordopt.o -L/usr/lib/R/lib -lR
installing to /home/syoh/R/x86_64-pc-linux-gnu-library/3.2/gconcordopt/libs
** R
** preparing package for lazy loading
** help
No man pages found in package  ‘gconcordopt’ 
*** installing help indices
** building package indices
** testing if installed package can be loaded
* DONE (gconcordopt)
syoh@skipper:~/tmp$ 

```

### Running ###


```
#!R

library(gconcordopt)
library(MASS)

p = 10
omega = diag(p)
omega[1,5] = omega[5,1] = .99
omega[2,6] = omega[6,2] = .99
sigma = solve(omega)

x = mvrnorm(200, mu=rep(0, p), Sigma=sigma)
omegahat = gconcordopt::concordista(x, lam=0.1)
eigen(omegahat)$values

```

with output


```
#!R

> library(gconcordopt)
> library(MASS)
> 
> p = 10
> omega = diag(p)
> omega[1,5] = omega[5,1] = .99
> omega[2,6] = omega[6,2] = .99
> sigma = solve(omega)
> 
> 
> x = mvrnorm(200, mu=rep(0, p), Sigma=sigma)
> omegahat = gconcordopt::concordista(x, lam=0.1)
> eigen(omegahat)$values
 [1] 1.78472360 1.76895084 1.08666961 0.98815492 0.98500556 0.97444126 0.97157266 0.95584726 0.01367573 0.01011393
> omega
      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
 [1,] 1.00 0.00    0    0 0.99 0.00    0    0    0     0
 [2,] 0.00 1.00    0    0 0.00 0.99    0    0    0     0
 [3,] 0.00 0.00    1    0 0.00 0.00    0    0    0     0
 [4,] 0.00 0.00    0    1 0.00 0.00    0    0    0     0
 [5,] 0.99 0.00    0    0 1.00 0.00    0    0    0     0
 [6,] 0.00 0.99    0    0 0.00 1.00    0    0    0     0
 [7,] 0.00 0.00    0    0 0.00 0.00    1    0    0     0
 [8,] 0.00 0.00    0    0 0.00 0.00    0    1    0     0
 [9,] 0.00 0.00    0    0 0.00 0.00    0    0    1     0
[10,] 0.00 0.00    0    0 0.00 0.00    0    0    0     1
> 
```