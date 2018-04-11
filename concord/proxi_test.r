library('R.matlab')

mat_data <- readMat('~/Dropbox/glasso/pMask.mat')
pMat <- mat_data$M
dimM <- dim(pMat)
print(dimM)

load('~/Dropbox/glasso/concord_results/RELATIONAL/lam7.3e-3_nz226.Rdata')
num_diag <- dim(inv_cov)[1]
idx <- which(inv_cov!=0,arr.ind=T)

num_nz <- dim(idx)[1]
num_cover <- 0.0
for(i in 1:num_nz){
	x<-idx[i,1]
	y<-idx[i,2]
	if(pMat[x,y]!=0){
		num_cover <- num_cover+1
	}
}
cover_rate <- num_cover/num_nz
print(cover_rate)

real_cover_rate <- (num_cover - num_diag)/(num_nz - num_diag)
print(real_cover_rate)

# load('lam5.5e-6_nz370_100.Rdata')
# a <- inv_cov
# load('lam5.5e-4_nz.Rdata')
# b <- inv_cov
# print(sum(a-b))
