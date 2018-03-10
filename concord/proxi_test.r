library('R.matlab')

mat_data <- readMat('~/Dropbox/glasso/pMask.mat')
pMat <- mat_data$M
dimM <- dim(pMat)
print(dimM)

load('~/Dropbox/glasso/concord_results/GAMBLING/lambda0.000055_nonzero370.Rdata')
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