library('R.matlab')
library('gconcord')
library(beepr)

mat_data <- readMat('~/Dropbox/glasso/data/ASD/JP_COR.mat')
data <- mat_data$X
dim <- dim(data)
print(dim)

inv_cov<-concord(data,0.8)
nonzero <- function(x) sum(x != 0)
print(nonzero(inv_cov)-dim[2])
writeMat('~/Dropbox/glasso/concord_results/ASD/0.8_.mat', M=inv_cov)
beep()

# mat_data <- readMat('~/Dropbox/glasso/pMask.mat')
# pMat <- mat_data$M

# inv_cov<-gconcordopt::concordista(data, lam=0.8, pMat=pMat)
# nonzero <- function(x) sum(x != 0)
# print(nonzero(inv_cov)-vec_dim)
# writeMat('~/Dropbox/glasso/concord_results/ASD/0.8_.mat', M=inv_cov)
# beep()